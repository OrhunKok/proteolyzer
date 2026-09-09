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
#
# The path is spelled out rather than reached through $HOME, because `docker exec
# -u` does not reliably set HOME -- and the earlier version of this, which used
# "$HOME/.ssh", wrote somewhere sshd never looks and failed with nothing but
# `Permission denied (publickey)`.
#
# chmod explicitly too. `mkdir -p` leaves an existing directory's mode alone, and
# the node image ships ~/.ssh as 755, so umask alone does not make it 700.
docker exec -i -u node "$container" sh -c '
    set -e
    mkdir -p /home/node/.ssh
    cat > /home/node/.ssh/authorized_keys
    chmod 700 /home/node/.ssh
    chmod 600 /home/node/.ssh/authorized_keys
' < "$key.pub"

# Read it back. A key that silently did not land is the failure this whole block
# exists to prevent, and it is invisible until ssh refuses.
if ! docker exec "$container" cat /home/node/.ssh/authorized_keys 2>/dev/null \
        | grep -qF "$(cut -d' ' -f2 < "$key.pub")"; then
    echo "cmux-attach: the public key is not in the container's authorized_keys." >&2
    echo "cmux-attach: container=$container key=$key.pub" >&2
    exit 1
fi

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

# Prove the login works before handing the connection to cmux. cmux reports a
# refused key as "the remote VM may have been paused, destroyed, or lost
# network", which is true of almost nothing and sends you looking in the wrong
# place; ssh's own message names the actual problem.
if ! ssh -o BatchMode=yes \
        -o IdentitiesOnly=yes \
        -o UserKnownHostsFile=/dev/null \
        -o StrictHostKeyChecking=no \
        -o ConnectTimeout=5 \
        -i "$key" -p "$port" node@127.0.0.1 true 2>/tmp/cmux-attach-ssh.$$; then
    echo "cmux-attach: ssh into the container failed. ssh said:" >&2
    sed 's/^/cmux-attach:   /' /tmp/cmux-attach-ssh.$$ >&2
    rm -f /tmp/cmux-attach-ssh.$$
    exit 1
fi
rm -f /tmp/cmux-attach-ssh.$$

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
#
# IdentitiesOnly=yes because `-i` only *adds* a key: ssh still offers everything
# in the agent first, and a well-stocked agent can exhaust MaxAuthTries before
# reaching the one key that would have worked. That failure also arrives as a
# bare `Permission denied (publickey)`.
exec cmux ssh "node@127.0.0.1" \
    --port "$port" \
    --identity "$key" \
    --name "$(basename "$repo")" \
    --no-forward-agent \
    --ssh-option IdentitiesOnly=yes \
    --ssh-option UserKnownHostsFile=/dev/null \
    --ssh-option StrictHostKeyChecking=no \
    --command "$remote_command"
