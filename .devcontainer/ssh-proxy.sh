#!/usr/bin/env bash
# Connect a socket to this project's container, for ssh's ProxyCommand.
#
#   ssh-proxy.sh <container-name> [port]
#
# Why this exists. Apple `container` does publish DNS records -- a container
# started with `container run --name dnstest` resolves immediately under the
# configured domain -- but a container created by `adevcontainer up` never gets
# one, verified on 2026-09-10 with the domain configured, the service restarted
# and the container recreated afterwards. So the name works for everything
# except the containers this repository actually makes.
#
# The name was wanted for one reason: every rebuild gets a fresh address, so a
# cmux workspace or a `cmux surface resume set` command pinned to an IP goes
# stale. A ProxyCommand solves that without DNS at all. ssh never resolves the
# alias -- it hands the hostname to this script, which looks the address up at
# connect time and pipes the socket through. The alias is stable; the address it
# reaches is whatever the container has right now.
set -euo pipefail

name="${1:?ssh-proxy.sh: need a container name}"
port="${2:-2222}"

# `container list` first: it is a table, so this depends on column order, but it
# costs nothing and does not need the container's shell. ID IMAGE OS ARCH STATE
# IP, so the address is the sixth field, with a prefix length to strip.
ip="$(container list 2>/dev/null \
    | awk -v n="$name" '$1 == n && $5 == "running" {print $6; exit}' \
    | cut -d/ -f1 || true)"

# Ask the container itself if that did not work -- slower, but it does not care
# what `container list` prints.
if [ -z "$ip" ]; then
    # `|| true` because pipefail plus set -e would otherwise abort here when
    # adevcontainer is missing, before the error message below can be printed --
    # which is exactly how this failed silently the first time.
    ip="$(adevcontainer exec --name "$name" -- hostname -i 2>/dev/null \
        | tr -d '\r' | awk '{print $1}' || true)"
fi

if [ -z "$ip" ]; then
    echo "ssh-proxy: no address for container '$name'. Is it running?" >&2
    echo "ssh-proxy:   ./.devcontainer/up.sh" >&2
    exit 1
fi

# `exec` so ssh owns the socket directly and there is no shell in the middle to
# mishandle EOF.
exec nc "$ip" "$port"
