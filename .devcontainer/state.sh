#!/usr/bin/env bash
# Move this container's state between machines, or between runtimes.
#
# The image definition travels in git and the workspace is a bind mount, so a
# new machine is `git clone`, `build.sh` and `up.sh`. What does not travel is the
# named volume, and it is the interesting part: /home/node/.state holds Claude
# Code's config, projects and memory under `claude/`, the GitHub login under
# `gh/`, and the shell history. Losing that on every machine change is what makes
# a portable environment feel unportable.
#
# Runs on the host.
#
#   ./.devcontainer/state.sh adopt                     take over the old VS Code
#                                                      container's state
#   ./.devcontainer/state.sh export                    without the GitHub token
#   ./.devcontainer/state.sh export --with-credentials  with it
#   ./.devcontainer/state.sh import [--with-credentials]
#   ./.devcontainer/state.sh migrate                   fold the old three into one
#   ./.devcontainer/state.sh export ~/somewhere.tar.gz
#
# `adopt` is the one for converting a project that already runs under VS Code and
# Docker. No arguments: it finds that project's old container by the label the
# devcontainer spec stamps on it, reads the Claude config, GitHub login and shell
# history straight out of it with `docker cp` -- so no volume names are involved,
# hashed or otherwise -- and writes them into this project's volume. The old
# container and its volumes are left untouched.
#
# Every export contains Claude Code's own credential either way; see the note in
# `export` below. Treat the archive as a secret. `adopt` writes no archive, which
# is the other reason to prefer it.
#
# The container need not be running -- but it must not be *running* for export or
# migrate, because a volume attaches to one container at a time.
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
        export|import|migrate|adopt) action="$1" ;;
        --with-credentials) with_credentials=1 ;;
        # Derived, not a line range: this printed `2,19p` while the header was
        # 19 lines, and grew silently wrong the moment the header did -- `--help`
        # cut off mid-sentence and dropped everything after it. Print the comment
        # block after the shebang and stop at the first line that is not one.
        -h|--help) awk 'NR==1{next} /^#/{sub(/^# ?/, ""); print; next} {exit}' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) archive="$1" ;;
    esac
    shift
done

if [ -z "$action" ]; then
    echo "state.sh: say 'adopt', 'export', 'import' or 'migrate'. --help for the rest." >&2
    exit 1
fi

# Runtimes keep separate volume stores, so this is also how you move state
# between them:
#
#   RUNTIME=docker    ./.devcontainer/state.sh export
#   RUNTIME=container ./.devcontainer/state.sh import
#
# Export reads and does not write, so the old setup stays intact until the new
# one is trusted.
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

# Volume mounts are not marked :ro anywhere below, and that is Apple
# `container`'s constraint rather than a preference. Mounting a named volume
# read-only to a second container is an open request -- apple/container#889 --
# and asking for it produces `VZErrorDomain Code=2 "The storage device
# attachment is invalid."`, which names the storage layer and not the flag. A
# host bind mount keeps :ro; that is a different mechanism and works.
#
# The same constraint means a volume already attached to a running container
# cannot attach to a second one, so `export` and `migrate` want the project's
# container stopped. run_or_hint says so when the attach fails.
run_or_hint() {
    if ! "$rt" "$@"; then
        echo >&2
        echo "state.sh: the runtime refused to attach a volume." >&2
        echo "state.sh: a volume attaches to one container at a time, so stop" >&2
        echo "state.sh: this project's container and try again:" >&2
        echo "state.sh:   adevcontainer stop" >&2
        exit 1
    fi
}

volume="claude-code-state-$base"

# The volume names from before everything moved under one. `migrate` reads these,
# and only matches a project whose old volumes were keyed on the folder basename.
# When they are not -- Anthropic's template keys them on ${devcontainerId}, so
# they come out as hashes -- use `adopt`, which asks the old container instead of
# guessing at names.
old_history="claude-code-bashhistory-$base"
old_config="claude-code-config-$base"
old_gh="claude-code-gh-$base"

case "$action" in
export)
    volume_exists "$volume" || {
        echo "state.sh: no volume $volume. Has the container ever started?" >&2
        exit 1
    }

    # `--with-credentials` governs the gh/ subtree and nothing else. This said
    # "an ordinary export carries no token", which is not true and was worth
    # correcting: claude/.credentials.json is Claude Code's own OAuth credential
    # and it is inside config/, so **every** export contains it.
    #
    # Left that way on purpose rather than fixed, because carrying the Claude
    # login is the point of an export when converting a project -- excluding it
    # would mean re-authenticating on the other side of every migration. The
    # consequence is the part that needs saying: any archive this writes is a
    # secret, `--with-credentials` or not.
    exclude=(--exclude=./gh)
    echo "state.sh: NOTE -- $archive will contain Claude Code's own credential"
    echo "state.sh: (claude/.credentials.json). Treat the file as a secret and"
    echo "state.sh: delete it once the other side is up. Do not leave it in a"
    echo "state.sh: synced folder -- the repository directory often is one."
    if [ "$with_credentials" -eq 1 ]; then
        exclude=()
        echo "state.sh: WARNING -- this also includes gh/, so $archive will hold a"
        echo "state.sh: GitHub token in plaintext. Move it as you would a key and"
        echo "state.sh: delete it after. \`gh auth login\` is one command."
    fi

    mkdir -p "$(dirname "$archive")"
    run_or_hint run --rm \
        -v "$volume:/v" \
        -v "$(dirname "$archive"):/out" \
        alpine tar czf "/out/$(basename "$archive")" \
            --numeric-owner ${exclude[@]+"${exclude[@]}"} -C /v .
    echo "state.sh: wrote $archive"
    ;;

import)
    [ -f "$archive" ] || { echo "state.sh: no such archive: $archive" >&2; exit 1; }
    ensure_volume "$volume"

    # An archive made without --with-credentials simply has no gh/ in it, so
    # there is nothing to decline here and nothing to half-restore.
    run_or_hint run --rm \
        -v "$volume:/v" \
        -v "$(dirname "$archive"):/in:ro" \
        alpine tar xzf "/in/$(basename "$archive")" --numeric-owner -C /v
    echo "state.sh: restored into $volume"
    echo "state.sh: the directory must stay named '$base' for the container to find it."
    ;;

migrate)
    # One-off, for a container that predates the single volume: fold
    # bashhistory, config and gh into .state/{history,claude,gh}. The old volumes
    # are mounted read-only and left in place, so this is repeatable and undone
    # by deleting the new volume.
    args=(--rm)
    found=0
    for pair in "$old_history:hist" "$old_config:config" "$old_gh:gh"; do
        name="${pair%%:*}"
        if volume_exists "$name"; then
            args+=(-v "$name:/old/${pair##*:}")
            found=1
            echo "state.sh: will read $name"
        else
            echo "state.sh: no $name; skipping."
        fi
    done
    [ "$found" -eq 1 ] || { echo "state.sh: nothing to migrate." >&2; exit 1; }

    ensure_volume "$volume"
    args+=(-v "$volume:/new")

    # 1000:1000 rather than a name: the alpine doing the copying has no `node`
    # user, and numeric ownership is what the container reads it back as.
    run_or_hint run "${args[@]}" alpine sh -c '
        set -e
        mkdir -p /new/claude /new/gh
        [ -d /old/config ] && cp -a /old/config/. /new/claude/ || true
        [ -d /old/gh ] && cp -a /old/gh/. /new/gh/ || true
        [ -f /old/hist/.bash_history ] && cp -a /old/hist/.bash_history /new/history || true
        chown -R 1000:1000 /new
    '
    echo "state.sh: folded into $volume; the old volumes are untouched."
    echo "state.sh: pick it up with: ./.devcontainer/build.sh && REBUILD=1 ./.devcontainer/up.sh"
    ;;

adopt)
    # Take over the state of this project's old VS Code / Docker devcontainer,
    # without being told anything about it.
    #
    # `docker cp` reads paths out of a container, volume-backed ones included, so
    # this needs no volume names -- which is the entire difficulty otherwise:
    # Anthropic's template keys its volumes on ${devcontainerId}, so they come out
    # as hashes that match nothing you could guess and nothing this script could
    # default to.
    #
    # The old container is found by the label the devcontainer spec already
    # stamps on it, so "which container" is not a question either. It is read
    # while stopped and never written to, and the old volumes are never touched:
    # the only thing this writes is the new volume, so a bad guess costs nothing
    # but a `container volume delete`.
    command -v docker >/dev/null 2>&1 || {
        echo "state.sh: adopt reads the old container through docker, which is not on PATH." >&2
        exit 1
    }
    docker info >/dev/null 2>&1 || { echo "state.sh: Docker is not running." >&2; exit 1; }

    old="${OLD_CONTAINER:-$(docker ps -a \
        --filter "label=devcontainer.local_folder=$repo" \
        --format '{{.ID}}' | head -1)}"
    if [ -z "$old" ]; then
        echo "state.sh: no Docker devcontainer is labelled with $repo." >&2
        echo "state.sh: Docker knows these, with the folder each was opened from:" >&2
        docker ps -a --format '  {{.ID}}  {{.Names}}  {{index .Labels "devcontainer.local_folder"}}' >&2
        echo "state.sh: if this project used to live under another path, name it:" >&2
        echo "state.sh:   OLD_CONTAINER=<id> $0 adopt" >&2
        exit 1
    fi
    echo "state.sh: adopting from container $old"

    staging="$(mktemp -d)"
    trap 'rm -rf "$staging"' EXIT
    found=0

    # Each of these is where the thing lives in a stock Claude devcontainer. A
    # miss is printed rather than swallowed: a silent partial adopt is the defect
    # worth avoiding -- the container comes up, and one of the three is missing.
    take() {
        if docker cp "$old:$1" "$staging/$2" >/dev/null 2>&1; then
            echo "state.sh: took $1"
            found=1
        else
            echo "state.sh: no $1 in that container; skipping"
        fi
    }
    take /home/node/.claude claude
    take /home/node/.config/gh gh
    # Two spellings, and the template's own is the second. Only the first that
    # exists is taken, so they cannot clobber each other.
    [ -e "$staging/history" ] || take /home/node/.bash_history history
    [ -e "$staging/history" ] || take /commandhistory/.bash_history history

    [ "$found" -eq 1 ] || {
        echo "state.sh: that container held none of the three; nothing to adopt." >&2
        exit 1
    }

    # Written with $rt, read with docker: that is the runtime crossing, and it is
    # why this is one command rather than an export and an import.
    ensure_volume "$volume"
    run_or_hint run --rm \
        -v "$volume:/new" \
        -v "$staging:/in" \
        alpine sh -c '
            set -e
            mkdir -p /new/claude /new/gh
            [ -d /in/claude ] && cp -a /in/claude/. /new/claude/ || true
            [ -d /in/gh ] && cp -a /in/gh/. /new/gh/ || true
            [ -f /in/history ] && cp -a /in/history /new/history || true
            chown -R 1000:1000 /new
        '
    echo "state.sh: adopted into $volume using $rt; the old container is untouched."
    echo "state.sh: start it with: ./.devcontainer/build.sh && ./.devcontainer/up.sh"
    ;;
esac
