#!/bin/bash
# Make the mounted volumes writable by node. Idempotent, run from
# postStartCommand ahead of everything else, and behind sudo for the same reason
# init-firewall.sh is.
#
# Docker seeds a named volume from the image path it covers, ownership included,
# so on Docker this finds everything already correct and does nothing. Apple
# `container` does not — apple/container#729, open since October 2025 — so all
# three volumes arrive empty and root-owned, `gh auth login` fails on a bare
# `Permission denied`, and Claude cannot write its own config directory. The
# Dockerfile creates these paths node-owned precisely so a seeding runtime
# inherits that; this is the same guarantee for a runtime that does not seed.
set -euo pipefail

for path in /commandhistory /home/node/.claude /home/node/.config/gh; do
    [ -d "$path" ] || continue
    owner="$(stat -c '%U' "$path")"
    if [ "$owner" != "node" ]; then
        echo "fix-volume-perms: $path was owned by $owner; giving it to node"
        chown -R node:node "$path"
    fi
done
