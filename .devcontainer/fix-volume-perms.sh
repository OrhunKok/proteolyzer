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

state=/home/node/.state

if [ -d "$state" ]; then
    owner="$(stat -c '%U' "$state")"
    if [ "$owner" != node ]; then
        echo "fix-volume-perms: $state was owned by $owner; giving it to node"
        chown -R node:node "$state"
    fi
    # The subdirectories exist in the image so a seeding runtime inherits them.
    # A runtime that does not seed gives an empty volume, so make them here too.
    mkdir -p "$state/claude/projects/-workspace/memory" "$state/gh"
    chown -R node:node "$state/claude" "$state/gh"
fi
