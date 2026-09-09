#!/usr/bin/env bash
# Open this repository's devcontainer as a cmux *remote workspace* rather than a
# pane with a shell in it.
#
# cmux has no devcontainer support and does not need any: `cmux ssh` is a
# first-class workspace type with a Linux-side daemon, and a container running
# sshd is a remote host like any other. cmux's own Docker integration suites
# attach to a container exactly this way -- ephemeral published port, throwaway
# key, host key checking off -- so this is a tested configuration rather than a
# clever one.
#
# What it buys over `up.sh` (which is `devcontainer exec`, and still the right
# tool for a quick shell):
#
#   - the `cmux` CLI works *inside* the container. The bootstrap installs a
#     wrapper at ~/.cmux/bin/cmux and pins CMUX_SOCKET_PATH to a reverse-
#     forwarded port, so `cmux notify` from a Claude Code hook reaches the app.
#   - PTY sessions survive cmux quitting and restarting, and reconnect.
#   - the sidebar gets real metadata, and dragging an image in uploads over sftp.
#   - browser panes egress from the container, so they are inside the firewall.
#
# This runs on the host, not inside the container.
#
#   ./.devcontainer/cmux-attach.sh                 workspace with a shell in /workspace
#   ./.devcontainer/cmux-attach.sh claude          workspace that starts Claude Code
#   REBUILD=1 ./.devcontainer/cmux-attach.sh       rebuild the container first
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
key="${CMUX_DEVCONTAINER_KEY:-$HOME/.ssh/cmux-devcontainer}"

for tool in cmux devcontainer docker; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "cmux-attach: $tool is not on PATH." >&2
        case "$tool" in
            cmux) echo "cmux-attach:   brew install --cask cmux" >&2 ;;
            devcontainer) echo "cmux-attach:   npm install -g @devcontainers/cli" >&2 ;;
            docker) echo "cmux-attach:   start Docker Desktop or OrbStack" >&2 ;;
        esac
        exit 1
    fi
done

# postStartCommand starts sshd, so this is also what guarantees it is running.
log="$(mktemp)"
trap 'rm -f "$log"' EXIT
# shellcheck disable=SC2086
devcontainer up --workspace-folder "$repo" ${REBUILD:+--remove-existing-container} | tee "$log"

container="$(sed -n 's/.*"containerId"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' "$log" | tail -1)"
if [ -z "$container" ]; then
    # The CLI labels every container it makes with the folder it came from.
    container="$(docker ps -q --filter "label=devcontainer.local_folder=$repo" | head -1)"
fi
rm -f "$log"
trap - EXIT   # the script ends in exec, which would never reach the trap

if [ -z "$container" ]; then
    echo "cmux-attach: could not work out which container is this workspace's." >&2
    exit 1
fi

# A key of its own, not the one that talks to GitHub: this authenticates a
# loopback hop into a container, and giving that a key with any other reach is
# how a sandbox stops being one.
if [ ! -f "$key" ]; then
    mkdir -p "$(dirname "$key")"
    ssh-keygen -t ed25519 -N '' -C cmux-devcontainer -f "$key"
fi

# Written on every attach rather than mounted: a mount that has to exist is a
# container that will not start on the day the file is missing.
docker exec -i -u node "$container" sh -c '
    umask 077
    mkdir -p "$HOME/.ssh"
    cat > "$HOME/.ssh/authorized_keys"
' < "$key.pub"

# Idempotent, and covers a container that was already up from before sshd was
# part of this image.
docker exec -u node "$container" sudo /usr/local/bin/start-sshd.sh

port="$(docker port "$container" 2222/tcp | head -1 | awk -F: '{print $NF}')"
if [ -z "$port" ]; then
    echo "cmux-attach: port 2222 is not published. Rebuild: REBUILD=1 $0" >&2
    exit 1
fi

if command -v nc >/dev/null 2>&1; then
    for _ in $(seq 1 30); do
        nc -z 127.0.0.1 "$port" 2>/dev/null && break
        sleep 0.2
    done
fi

remote_command="cd /workspace"
if [ "$#" -gt 0 ]; then
    remote_command="cd /workspace && $*"
fi

# Host key checking off, deliberately: start-sshd.sh generates a host key per
# container and the port is a fresh ephemeral one each rebuild, so known_hosts
# would do nothing but reject a container it has seen before under a port some
# other container used. What bounds this is the port being on loopback and the
# key being one this script made.
#
# --no-forward-agent for the reason the container holds its own gh credential:
# it is a sandbox, and the host's ssh agent is not part of what it gets.
exec cmux ssh "node@127.0.0.1" \
    --port "$port" \
    --identity "$key" \
    --name "$(basename "$repo")" \
    --no-forward-agent \
    --ssh-option UserKnownHostsFile=/dev/null \
    --ssh-option StrictHostKeyChecking=no \
    --command "$remote_command"
