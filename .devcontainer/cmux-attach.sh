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
            # cmux puts its CLI on the PATH of terminals *it* spawns, not on the
            # system one, so "not on PATH" usually means a plain Terminal.app
            # window rather than a missing install. Suggesting `brew install`
            # first sends you to reinstall something you already have.
            cmux)
                echo "cmux-attach:   run this from a cmux tab -- cmux only puts" >&2
                echo "cmux-attach:   its CLI on the PATH of terminals it starts." >&2
                echo "cmux-attach:   Not installed at all? brew install --cask cmux" >&2
                ;;
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

# A stable name for the container, which is the whole point: every rebuild gets a
# fresh address, so a cmux workspace or a `cmux surface resume set` command
# pinned to an IP goes stale the next time you rebuild.
#
# Not via DNS. Apple `container` does publish records -- `container run --name
# dnstest` resolved immediately under the configured domain -- but a container
# created by `adevcontainer up` never gets one, checked on 2026-09-10 with the
# domain in `container system property list`, the service restarted and the
# container recreated after it. So the name works for everything except the
# containers this repository makes.
#
# An ssh alias does not need DNS. ssh matches the alias literally, never resolves
# it, and hands the hostname to ssh-proxy.sh, which looks the address up at
# connect time. The alias is stable; the address behind it is whatever the
# container has right now. That is strictly better than a DNS record, which would
# still be one rebuild behind between restarts.
host="${CMUX_DEVCONTAINER_HOST:-$(basename "$repo").${CONTAINER_DNS_DOMAIN:-adevcontainers.local}}"
ssh_config="$HOME/.ssh/config"
marker_begin="# BEGIN cmux-devcontainer $host"
marker_end="# END cmux-devcontainer $host"

# Prepended, not appended, and that matters: ssh takes the *first* value it sees
# for each keyword, so a `Host *` block earlier in the file would win on
# IdentityFile and the right key would never be offered.
mkdir -p "$HOME/.ssh"
touch "$ssh_config"
if ! grep -qF "$marker_begin" "$ssh_config"; then
    block="$(mktemp)"
    {
        echo "$marker_begin"
        echo "# Written by .devcontainer/cmux-attach.sh. Delete this block to opt out."
        echo "Host $host"
        echo "    User node"
        echo "    ProxyCommand \"$repo/.devcontainer/ssh-proxy.sh\" \"$(basename "$repo")\" 2222"
        echo "    IdentityFile $key"
        echo "    IdentitiesOnly yes"
        echo "    StrictHostKeyChecking no"
        echo "    UserKnownHostsFile /dev/null"
        echo "$marker_end"
        echo
        cat "$ssh_config"
    } > "$block"
    mv "$block" "$ssh_config"
    chmod 600 "$ssh_config"
    echo "cmux-attach: added an ssh alias for $host to $ssh_config"
fi

echo "cmux-attach: using $host"

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
