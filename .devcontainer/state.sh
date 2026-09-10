#!/usr/bin/env bash
# Move this container's state between machines.
#
# The image definition travels in git and the workspace is a bind mount, so a
# new machine is `git clone` and `up.sh`. What does not travel is the three
# named volumes, and one of them is the interesting one: /home/node/.claude
# holds settings, project state and the agent's memory directory. Losing that
# on every machine change is the part that makes a portable environment feel
# unportable.
#
# Runs on the host.
#
#   ./.devcontainer/state.sh export                  everything but the credential
#   ./.devcontainer/state.sh export --with-credentials
#   ./.devcontainer/state.sh import [--with-credentials]
#   ./.devcontainer/state.sh export ~/somewhere.tar.gz
#
# The container does not need to be running. `docker volume` is enough.
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
base="$(basename "$repo")"
archive="$repo/devcontainer-state.tar.gz"
with_credentials=0
action=""

while [ "$#" -gt 0 ]; do
    case "$1" in
        export|import) action="$1" ;;
        --with-credentials) with_credentials=1 ;;
        -h|--help) sed -n '2,18p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) archive="$1" ;;
    esac
    shift
done

if [ -z "$action" ]; then
    echo "state.sh: say 'export' or 'import'. --help for the rest." >&2
    exit 1
fi

# Runtimes keep separate volume stores, so this is also how you move state
# between them -- which is the migration you want the day you switch:
#
#   RUNTIME=docker    ./.devcontainer/state.sh export
#   RUNTIME=container ./.devcontainer/state.sh import
#
# The Docker volumes are read, not touched, so the old setup stays intact until
# you are satisfied the new one works.
rt="${RUNTIME:-}"
if [ -z "$rt" ]; then
    if command -v container >/dev/null 2>&1; then
        rt=container
    elif command -v docker >/dev/null 2>&1; then
        rt=docker
    else
        echo "state.sh: no container runtime found." >&2
        exit 1
    fi
fi

if ! command -v "$rt" >/dev/null 2>&1; then
    echo "state.sh: RUNTIME=$rt is not on PATH." >&2
    exit 1
fi

if [ "$rt" = docker ] && ! docker info >/dev/null 2>&1; then
    echo "state.sh: Docker is not running." >&2
    exit 1
fi

echo "state.sh: using $rt"

# `container volume` has create/delete/prune/list but no `inspect`, and its
# `create` errors on an existing volume where Docker's is idempotent. `volume
# list --quiet` is spelled the same on both, so existence goes through that.
volume_exists() {
    "$rt" volume list --quiet 2>/dev/null | grep -qx "$1"
}

ensure_volume() {
    volume_exists "$1" || "$rt" volume create "$1" >/dev/null
}

# Volume names are keyed on the directory basename, the same as devcontainer.json
# does it -- so the directory has to be named the same on the far machine for
# these to land where the container will look for them.
history_volume="claude-code-bashhistory-$base"
config_volume="claude-code-config-$base"
gh_volume="claude-code-gh-$base"

case "$action" in
export)
    args=(--rm)
    for pair in "$history_volume:bashhistory" "$config_volume:config"; do
        volume="${pair%%:*}"
        if volume_exists "$volume"; then
            args+=(-v "$volume:/v/${pair##*:}:ro")
        else
            echo "state.sh: no volume $volume yet; skipping." >&2
        fi
    done

    if [ "$with_credentials" -eq 1 ]; then
        if volume_exists "$gh_volume"; then
            args+=(-v "$gh_volume:/v/gh:ro")
            echo "state.sh: WARNING -- including $gh_volume puts a GitHub token"
            echo "state.sh: in $archive as plaintext. Move it as you would a key,"
            echo "state.sh: and delete it after. \`gh auth login\` is one command."
        else
            echo "state.sh: no volume $gh_volume yet; skipping." >&2
        fi
    fi

    mkdir -p "$(dirname "$archive")"
    args+=(-v "$(dirname "$archive"):/out")
    "$rt" run "${args[@]}" alpine \
        tar czf "/out/$(basename "$archive")" --numeric-owner -C /v .
    echo "state.sh: wrote $archive"
    ;;

import)
    [ -f "$archive" ] || { echo "state.sh: no such archive: $archive" >&2; exit 1; }

    args=(--rm)
    ensure_volume "$history_volume"
    ensure_volume "$config_volume"
    args+=(-v "$history_volume:/v/bashhistory" -v "$config_volume:/v/config")

    if [ "$with_credentials" -eq 1 ]; then
        ensure_volume "$gh_volume"
        args+=(-v "$gh_volume:/v/gh")
    fi

    args+=(-v "$(dirname "$archive"):/in:ro")
    # Anything in the archive without a volume mounted over its path lands in the
    # throwaway container layer and goes away with it -- which is how a `gh`
    # directory in the tarball is declined rather than half-restored.
    "$rt" run "${args[@]}" alpine \
        tar xzf "/in/$(basename "$archive")" --numeric-owner -C /v
    echo "state.sh: restored into $history_volume, $config_volume$([ "$with_credentials" -eq 1 ] && echo ", $gh_volume")"
    echo "state.sh: the directory must stay named '$base' for the container to find these."
    ;;
esac
