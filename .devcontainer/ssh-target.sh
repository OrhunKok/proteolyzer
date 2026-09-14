# Bring this repository's devcontainer up and make it reachable over ssh, then
# say where it is. Sourced, not run -- it sets shell variables for its caller:
#
#   repo       absolute path of the repository on the Mac
#   container  the managed container's name, from the runtime rather than guessed
#   key        the ssh identity that authenticates into it
#   target     the DNS name if this Mac resolves it, else the address
#
# This is everything a frontend needs and nothing any particular frontend wants,
# which is the reason it is a file of its own: `cmux-attach.sh` and
# `orca-target.sh` both need all of it, and the alternative was two copies of a
# hundred lines whose every paragraph records a bug that already happened once.
# This repository has been through that with `CoordinatesMapping` and does not
# need to do it again.
#
# Apple `container` only. Every container gets its own address reachable from the
# host, so there is no published port, no `docker port` and no loopback juggling
# -- which is most of why this is worth doing on that runtime. On Docker use
# `up.sh`; the config publishes no port there.
#
# Callers are expected to have `set -euo pipefail` already; this does not set it,
# because a sourced file changing its caller's shell options is a surprise.

# Message prefix. `$0` is still the *calling* script inside a sourced file, so
# this reads `cmux-attach:` or `orca-target:` without either having to say which.
_st_me="${0##*/}"
_st_me="${_st_me%.sh}"

# These drive macOS-side tooling -- `container`, `adevcontainer`, whichever
# frontend called -- so running one *inside* the container is a mistake worth
# naming. Left uncaught the symptom is "not on PATH" plus an invitation to
# `brew install` it, on Linux, which sends you somewhere with no exit.
if [ "$(uname -s)" = Linux ]; then
    printf '%s\n' \
        "$_st_me: this runs on the Mac, not inside the container." \
        "$_st_me: \`exit\` back to the host first, or open a terminal there." >&2
    exit 1
fi

# The frontend's own tool check, if it has one, runs here: after the guard above
# and before anything is built.
#
# Both halves of that matter. A frontend checked *before* the guard turns "you
# are inside the container" into "cmux is not on PATH, try brew install" on
# Linux, which is the trap the guard exists for -- `cmux` is genuinely absent in
# a container that has not been attached to yet, so it is the check most likely
# to fire there and mislead. And a frontend checked *after* the build discovers
# it is missing two minutes into an `adevcontainer up`.
#
# A function rather than a list of names because the useful part of these checks
# is the hint, and cmux's is four lines about how its CLI reaches a terminal.
if declare -F st_frontend_check >/dev/null 2>&1; then
    st_frontend_check
fi

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# One identity for every project, and a dedicated one: this authenticates a hop
# into a sandbox, and giving that a key with any other reach -- the GitHub key,
# say -- is how a sandbox stops being one.
key="${CMUX_DEVCONTAINER_KEY:-$HOME/.ssh/cmux-devcontainer}"

for tool in adevcontainer container; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "$_st_me: $tool is not on PATH." >&2
        case "$tool" in
            adevcontainer) echo "$_st_me:   brew install wcgomes/tap/adevcontainer" >&2 ;;
            container) echo "$_st_me:   see github.com/apple/container (macOS 26+)" >&2 ;;
        esac
        echo "$_st_me: on Docker, use ./.devcontainer/up.sh instead." >&2
        exit 1
    fi
done

cd "$repo"

# Every `adevcontainer exec` below is given --name, and that is not tidiness.
# With more than one managed container running, `exec` without it opens an
# interactive "Select a container:" picker -- so a script that does four execs
# stops four times, and whichever entry happens to be highlighted is the
# container it acts on. That is how this ended up writing a key into one
# project, reading sshd's pid from a second, and reporting a third one's address
# as this project's: the alias for `proteolyzer` pointed at `pinpoint`.
#
# The name comes from `containerId:` in the tool's own output rather than being
# guessed from the directory, so it stays right even if `name` in
# devcontainer.json and the folder ever disagree.
_st_log="$(mktemp)"
trap 'rm -f "$_st_log"' EXIT
if [ -n "${REBUILD:-}" ]; then
    # `rebuild --name` selects an *existing* container, so the name has to be the
    # one that exists -- not the one devcontainer.json would create. Those differ
    # whenever `name` has been edited or the folder renamed, and the failure is
    # `No managed container named ...` while the container is sitting right there
    # under another name. `list` knows which container belongs to this folder;
    # ask it. The leading slash in the match keeps `notpinpoint` from answering
    # for `pinpoint`.
    _st_existing="$(adevcontainer list 2>/dev/null \
        | awk -v b="/$(basename "$repo")" \
            'NR>1 { n=length(b); if (substr($NF, length($NF)-n+1) == b) { print $1; exit } }' \
        || true)"

    if [ -n "$_st_existing" ]; then
        adevcontainer rebuild --name "$_st_existing" 2>&1 | tee "$_st_log"
    else
        echo "$_st_me: no container for this folder yet; creating one." >&2
        adevcontainer up 2>&1 | tee "$_st_log"
    fi
else
    adevcontainer up 2>&1 | tee "$_st_log"
fi

container="$(sed -n 's/.*containerId:[[:space:]]*\([A-Za-z0-9_.-][A-Za-z0-9_.-]*\).*/\1/p' "$_st_log" | tail -1)"
rm -f "$_st_log"
trap - EXIT
if [ -z "$container" ]; then
    container="$(basename "$repo")"
    echo "$_st_me: no containerId in the output; assuming '$container'." >&2
fi

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
adevcontainer exec --name "$container" -- sh -c '
    set -e
    mkdir -p /home/node/.ssh
    printf "%s\n" "$1" > /home/node/.ssh/authorized_keys
    chmod 700 /home/node/.ssh
    chmod 600 /home/node/.ssh/authorized_keys
' sh "$(cat "$key.pub")"

# Read it back. A key that silently did not land is invisible until ssh refuses.
if ! adevcontainer exec --name "$container" -- cat /home/node/.ssh/authorized_keys 2>/dev/null \
        | grep -qF "$(cut -d' ' -f2 < "$key.pub")"; then
    echo "$_st_me: the public key is not in the container's authorized_keys." >&2
    exit 1
fi

adevcontainer exec --name "$container" -- sudo /usr/local/bin/start-sshd.sh

# Built from the container's name rather than the directory's. They agree when
# `name` in devcontainer.json matches the folder, and when they do not, the alias
# should follow the thing it actually reaches.
host="${CMUX_DEVCONTAINER_HOST:-$container.${CONTAINER_DNS_DOMAIN:-adevcontainers.local}}"

# The address, asked of the container rather than read out of a CLI table.
#
# This used to parse `container list` on column position -- ID IMAGE OS ARCH
# STATE IP, so field six, minus a prefix length. That worked and was a hostage
# to a format nobody promised to keep. `hostname -i` inside the container is
# authoritative, costs one exec on a path that runs once per attach, and cannot
# be broken by a column being added.
ip="$(adevcontainer exec --name "$container" -- hostname -i 2>/dev/null | tr -d '\r' | awk '{print $1}' || true)"
if [ -z "$ip" ]; then
    echo "$_st_me: could not read the container's address." >&2
    echo "$_st_me: it should be running by now -- check \`container list\`." >&2
    exit 1
fi

# Which of the two to hand the frontend. Nothing is written to ~/.ssh/config: a
# file every project appends a block to is shared state, and one container has no
# business knowing which others exist. Whatever a frontend needs, it is passed.
#
# The name is preferred when the Mac can actually resolve it, and the address is
# used when it cannot, so a machine without the DNS domain set up is not stuck.
# `dscacheutil` rather than `dig @127.0.0.1 -p 2053`: the second proves the
# record exists, the first proves *this Mac* will find it, and only the second
# question decides whether an attach works.
target="$ip"
if dscacheutil -q host -a name "$host" 2>/dev/null | grep -q '^ip_address:'; then
    target="$host"
fi
echo "$_st_me: $container at $target"

# Prove the login before handing the connection to anything else. cmux reports a
# refused key as "the remote VM may have been paused, destroyed, or lost network"
# -- true of almost nothing, and it sends you looking in the wrong place. Orca is
# clearer about it but still after the fact.
#
# IdentitiesOnly because `-i` only *adds* a key: ssh offers the agent's first,
# and a full agent can exhaust MaxAuthTries before reaching the one that works.
_st_err="$(mktemp)"
if ! ssh -o BatchMode=yes \
        -o IdentitiesOnly=yes \
        -o UserKnownHostsFile=/dev/null \
        -o StrictHostKeyChecking=no \
        -o ConnectTimeout=10 \
        -i "$key" -p 2222 "node@$target" true 2>"$_st_err"; then
    echo "$_st_me: ssh into the container failed. ssh said:" >&2
    sed "s/^/$_st_me:   /" "$_st_err" >&2
    rm -f "$_st_err"
    exit 1
fi
rm -f "$_st_err"
