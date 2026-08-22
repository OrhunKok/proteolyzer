#!/usr/bin/env bash
# Find claude-shared and run its installer. Copy this into a project's
# .devcontainer/ and call it from postCreateCommand:
#
#   "postCreateCommand": "bash .devcontainer/link-ecosystem.sh"
#
# This is the one piece that cannot live in claude-shared itself: something has
# to locate it before anything in it can run.
#
# Ordered by preference, and every tier is optional -- a project with none of
# them still builds, it is just not linked, and says so.

set -uo pipefail

REPO_URL="${CLAUDE_SHARED_URL:-https://github.com/OrhunKok/claude-shared}"
CACHE="$HOME/.claude-shared"

for candidate in \
    /workspace-shared \
    /workspace/claude-shared \
    "$CACHE"
do
    if [ -f "$candidate/install.sh" ]; then
        echo "link-ecosystem: using $candidate"
        exec bash "$candidate/install.sh"
    fi
done

# Nothing local. Cloning needs a credential for a private repository, which may
# not be available this early, so this is a best effort rather than a
# requirement.
echo "link-ecosystem: no local claude-shared; trying to clone it."
if git clone -q "$REPO_URL" "$CACHE" 2>/dev/null && [ -f "$CACHE/install.sh" ]; then
    exec bash "$CACHE/install.sh"
fi

cat <<'EOF'
link-ecosystem: could not find or clone claude-shared, so this container is not
linked to the other projects. Nothing else is affected.

To fix it, either clone claude-shared next to this repository on the host, so
the devcontainer mount resolves:

    git clone https://github.com/OrhunKok/claude-shared

or clone it inside the container and re-run:

    git clone https://github.com/OrhunKok/claude-shared ~/.claude-shared
    bash ~/.claude-shared/install.sh
EOF
