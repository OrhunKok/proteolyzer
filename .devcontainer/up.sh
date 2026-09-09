#!/usr/bin/env bash
# Bring this repository's devcontainer up and open a shell in it -- the job the
# VS Code Dev Containers extension used to do, minus the editor.
#
# This runs on the host, not inside the container. Call it by hand, or from a
# cmux custom command or resume command; see README.md.
#
#   ./.devcontainer/up.sh             a zsh inside the container
#   ./.devcontainer/up.sh claude      straight into Claude Code
#   REBUILD=1 ./.devcontainer/up.sh   discard the existing container first
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if ! command -v devcontainer >/dev/null 2>&1; then
    echo "up.sh: the devcontainer CLI is not installed." >&2
    echo "up.sh:   npm install -g @devcontainers/cli" >&2
    exit 1
fi

if ! docker info >/dev/null 2>&1; then
    echo "up.sh: Docker is not running." >&2
    exit 1
fi

# Idempotent: reuses a container that is already up, and runs init-firewall.sh
# through postStartCommand when it had to start one. Output is left visible
# because a first build takes minutes and silence there reads as a hang.
# shellcheck disable=SC2086
devcontainer up --workspace-folder "$repo" ${REBUILD:+--remove-existing-container}

# `zsh`, not `zsh -l`: it matches the profile VS Code's terminal used, and it is
# the interactive shell that sources ~/.zshrc, which is where gh-auth.sh lives.
exec devcontainer exec --workspace-folder "$repo" "${@:-zsh}"
