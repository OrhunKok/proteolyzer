#!/usr/bin/env bash
# Open this repository's devcontainer as a cmux *remote workspace* rather than a
# pane with a shell in it.
#
# cmux has no devcontainer support and does not need any: `cmux ssh` is a
# first-class workspace type with a Linux-side daemon, and a container running
# sshd is a remote host like any other.
#
# What it buys over `up.sh`, which is still the right tool for a quick shell:
#
#   - the `cmux` CLI works *inside* the container, so `cmux notify` and
#     `cmux workspace status set` from a Claude Code hook reach the app
#   - terminals survive cmux quitting, and reconnect on relaunch
#   - sidebar metadata, and dragging a file into the pane uploads over sftp
#
# Apple `container` only. Every container gets its own address reachable from the
# host, so there is no published port, no `docker port` and no loopback juggling
# -- which is most of why this is worth doing on that runtime. On Docker use
# `up.sh`; the config publishes no port there.
#
#   ./.devcontainer/cmux-attach.sh             a workspace with a shell
#   ./.devcontainer/cmux-attach.sh claude      a workspace running Claude Code
#   REBUILD=1 ./.devcontainer/cmux-attach.sh   rebuild first
set -euo pipefail

# These drive macOS-side tooling -- cmux, `container`, `adevcontainer`, Docker
# on the Mac -- so running one *inside* the container is a mistake worth naming.
# Left uncaught the symptom is "cmux is not on PATH" plus an invitation to
# `brew install` it, on Linux, which sends you somewhere with no exit.
if [ "$(uname -s)" = Linux ]; then
    printf '%s\n' \
        "${0##*/}: this runs on the Mac, not inside the container." \
        "${0##*/}: \`exit\` back to the host first, or use another cmux tab." >&2
    exit 1
fi

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
key="${CMUX_DEVCONTAINER_KEY:-$HOME/.ssh/cmux-devcontainer}"

for tool in cmux adevcontainer container; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "cmux-attach: $tool is not on PATH." >&2
        case "$tool" in
            cmux) echo "cmux-attach:   brew install --cask cmux" >&2 ;;
            adevcontainer) echo "cmux-attach:   brew install wcgomes/tap/adevcontainer" >&2 ;;
            container) echo "cmux-attach:   see github.com/apple/container (macOS 26+)" >&2 ;;
        esac
        echo "cmux-attach: on Docker, use ./.devcontainer/up.sh instead." >&2
        exit 1
    fi
done

cd "$repo"
if [ -n "${REBUILD:-}" ]; then
    adevcontainer rebuild
else
    adevcontainer up
fi

# A key of its own, not the one that talks to GitHub: this authenticates a hop
# into a sandbox, and giving that a key with any other reach is how a sandbox
# stops being one.
if [ ! -f "$key" ]; then
    mkdir -p "$(dirname "$key")"
    ssh-keygen -t ed25519 -N '' -C cmux-devcontainer -f "$key"
fi

# The key goes in as an argument, not on stdin. `adevcontainer exec` without -i
# attaches no stdin, so a `cat >` inside reads EOF immediately and writes an
# empty authorized_keys -- which is exactly what the read-back below caught. An
# argument needs nothing of the CLI beyond running the command.
#
# Spelled-out paths, and both modes set explicitly: `mkdir -p` leaves an existing
# directory's mode alone and the node image ships ~/.ssh as 755, so a umask alone
# never makes it 700.
adevcontainer exec -- sh -c '
    set -e
    mkdir -p /home/node/.ssh
    printf "%s\n" "$1" > /home/node/.ssh/authorized_keys
    chmod 700 /home/node/.ssh
    chmod 600 /home/node/.ssh/authorized_keys
' sh "$(cat "$key.pub")"

# Read it back. A key that silently did not land is invisible until ssh refuses.
if ! adevcontainer exec -- cat /home/node/.ssh/authorized_keys 2>/dev/null \
        | grep -qF "$(cut -d' ' -f2 < "$key.pub")"; then
    echo "cmux-attach: the public key is not in the container's authorized_keys." >&2
    exit 1
fi

adevcontainer exec -- sudo /usr/local/bin/start-sshd.sh

# Where to reach the container. A name is preferred over an address for a reason
# that matters here: every rebuild gets a fresh IP, so a workspace or a
# `cmux surface resume set` command pinned to an address goes stale the next time
# you rebuild, while a name does not.
#
# Names come from Apple `container`'s embedded DNS, which needs two one-off steps
# on the Mac -- see README.md. The suffix is `adev.containers`, so this project
# is `proteolyzer.adev.containers`: it says what answers the query and what kind
# of thing answered. `.containers` is not a delegated TLD, so nothing on the
# public internet can shadow it or be shadowed by it -- which a suffix ending in
# a real TLD like `.net` could not promise.
#
# Unconfigured, this falls back to the address and everything still works.
host="${CMUX_DEVCONTAINER_HOST:-}"
if [ -z "$host" ]; then
    candidate="$(basename "$repo").${CONTAINER_DNS_DOMAIN:-adev.containers}"
    # dscacheutil rather than dig: it goes through macOS's resolver, which is
    # what /etc/resolver configures and therefore what actually decides.
    if dscacheutil -q host -a name "$candidate" 2>/dev/null | grep -q '^ip_address:'; then
        host="$candidate"
        echo "cmux-attach: using $host"
    fi
fi

if [ -z "$host" ]; then
    host="$(adevcontainer exec -- hostname -i | tr -d '\r' | awk '{print $1}')"
    if [ -z "$host" ]; then
        echo "cmux-attach: could not reach the container by name or address." >&2
        exit 1
    fi
    echo "cmux-attach: using $host (no DNS; see README.md to get a name)"
fi

# Prove the login before handing the connection to cmux, which reports a refused
# key as "the remote VM may have been paused, destroyed, or lost network" -- true
# of almost nothing, and it sends you looking in the wrong place.
#
# IdentitiesOnly because `-i` only *adds* a key: ssh offers the agent's first,
# and a full agent can exhaust MaxAuthTries before reaching the one that works.
if ! ssh -o BatchMode=yes \
        -o IdentitiesOnly=yes \
        -o UserKnownHostsFile=/dev/null \
        -o StrictHostKeyChecking=no \
        -o ConnectTimeout=10 \
        -i "$key" -p 2222 "node@$host" true 2>/tmp/cmux-attach-ssh.$$; then
    echo "cmux-attach: ssh into the container failed. ssh said:" >&2
    sed 's/^/cmux-attach:   /' /tmp/cmux-attach-ssh.$$ >&2
    rm -f /tmp/cmux-attach-ssh.$$
    exit 1
fi
rm -f /tmp/cmux-attach-ssh.$$

# sshd builds its own login environment, so what is on PATH there depends on
# image plumbing this script has no business assuming. npm's global bin -- where
# `claude` is -- is added here, single-quoted so $PATH expands in the container
# and not on the Mac. Belongs in the image too, and is, but a script that only
# works against a freshly built image is a script that fails at the worst time.
prelude='export PATH="$PATH:/usr/local/share/npm-global/bin"'

remote_command="$prelude && cd /workspace"
if [ "$#" -gt 0 ]; then
    remote_command="$prelude && cd /workspace && $*"
fi

# Host key checking off: start-sshd.sh generates a host key per container, so
# known_hosts could only ever reject a rebuild of the same workspace. What bounds
# this is the address being the runtime's own and the key being one this made.
#
# --no-forward-agent for the reason the container holds its own gh credential:
# it is a sandbox, and the host's ssh agent is not part of what it gets.
exec cmux ssh "node@$host" \
    --port 2222 \
    --identity "$key" \
    --name "$(basename "$repo")" \
    --no-forward-agent \
    --ssh-option IdentitiesOnly=yes \
    --ssh-option UserKnownHostsFile=/dev/null \
    --ssh-option StrictHostKeyChecking=no \
    --command "$remote_command"
