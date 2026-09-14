#!/usr/bin/env bash
# Move this container's state between machines, or between runtimes.
#
# The image definition travels in git and the workspace is a bind mount, so a
# new machine is `git clone`, `build.sh` and `up.sh`. What does not travel is the
# named volume, and it is the interesting part: /home/node/.state holds Claude
# Code's config, projects and memory under `claude/`, the GitHub login under
# `gh/`, the git identity under `git/`, and the shell history. Losing that on
# every machine change is what makes a portable environment feel unportable.
#
# `export` tars the volume whole and excludes rather than enumerates, so a new
# subdirectory on it travels without this script being told about it. `git/` is
# name and email, which are not secrets; `gh/` holds a token and is the one
# subtree left out unless asked for.
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

# These drive macOS-side tooling -- `container`, `adevcontainer`, Docker on the
# Mac -- so running one *inside* the container is a mistake worth naming. Left
# uncaught the symptom is a tool "not on PATH" plus an invitation to
# `brew install` it, on Linux, which sends you somewhere with no exit.
if [ "$(uname -s)" = Linux ]; then
    printf '%s\n' \
        "${0##*/}: this runs on the Mac, not inside the container." \
        "${0##*/}: \`exit\` back to the host first, or open a terminal there." >&2
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

adopt)
    # Bring this project's pre-existing Claude state into its volume, from
    # wherever it happens to live, and never clobber what is already there.
    #
    # This exists because copying .devcontainer into a project gave it an empty
    # volume and no route to its own history -- a repository used with Claude for
    # months opened with no sessions, and the only remedy was a hand-run sequence
    # per project. `up.sh` calls this once per volume.
    #
    # The target volume is never mounted here, and that is the whole design. It
    # is attached to the running container, and a volume attaches to one
    # container at a time (apple/container#889) -- so the obvious implementation,
    # mounting source and target into a throwaway container, fails with
    # `VZErrorDomain Code=2 "The storage device attachment is invalid."` every
    # single time up.sh calls it, because up.sh starts the container first. The
    # workspace is a bind mount, so the source goes out as a tarball there and is
    # unpacked by the container that already holds the target.
    ensure_volume "$volume"

    # Which container holds this folder. `list` is the authority; `name` in
    # devcontainer.json may have been edited and the folder may have been renamed.
    into="$(adevcontainer list 2>/dev/null \
        | awk -v b="/$(basename "$repo")" \
            'NR>1 { n=length(b); if (substr($NF, length($NF)-n+1) == b) { print $1; exit } }' \
        || true)"
    if [ -z "$into" ]; then
        echo "state.sh: no running container for $base; start it first:" >&2
        echo "state.sh:   ./.devcontainer/up.sh" >&2
        exit 1
    fi

    src_rt=""
    if volume_exists "$old_config" || volume_exists "$old_history"; then
        src_rt="$rt"
    elif [ "$rt" != docker ] && command -v docker >/dev/null 2>&1 \
         && docker info >/dev/null 2>&1 \
         && docker volume list --quiet 2>/dev/null | grep -qx "$old_config"; then
        src_rt=docker
    fi

    # Whether Docker could be asked at all decides whether "found nothing" means
    # anything. If it is not running the legacy volumes are simply invisible, and
    # marking the volume adopted on that basis retires the search forever, on
    # exactly the run least able to answer the question.
    docker_checked=1
    if [ "$rt" != docker ]; then
        if ! command -v docker >/dev/null 2>&1 || ! docker info >/dev/null 2>&1; then
            docker_checked=0
        fi
    fi

    if [ -z "$src_rt" ]; then
        if [ "$docker_checked" -eq 0 ]; then
            echo "state.sh: Docker is not running, so earlier state kept in its volumes" >&2
            echo "state.sh: could not be looked for. Not marking this done -- start Docker" >&2
            echo "state.sh: Desktop and run up.sh again if this project had history." >&2
            exit 0
        fi
        echo "state.sh: nothing to adopt for $base."
        adevcontainer exec --name "$into" -- touch /home/node/.state/.adopted >/dev/null 2>&1 || true
        exit 0
    fi

    echo "state.sh: adopting $base's earlier state from $src_rt"

    # Staged inside .devcontainer because that is under the bind mount, so the
    # container sees it at /workspace/.devcontainer without anything being
    # attached twice.
    stage_host="$repo/.devcontainer/.adopt.tgz"
    stage_guest="/workspace/.devcontainer/.adopt.tgz"
    rm -f "$stage_host"

    args=(--rm)
    if [ "$src_rt" = docker ]; then
        docker volume list --quiet | grep -qx "$old_history" && args+=(-v "$old_history:/old/hist:ro")
        docker volume list --quiet | grep -qx "$old_config" && args+=(-v "$old_config:/old/config:ro")
        docker volume list --quiet | grep -qx "$old_gh" && args+=(-v "$old_gh:/old/gh:ro")
        args+=(-v "$repo/.devcontainer:/out")
        docker run "${args[@]}" alpine tar czf /out/.adopt.tgz --numeric-owner -C /old .
    else
        volume_exists "$old_history" && args+=(-v "$old_history:/old/hist")
        volume_exists "$old_config" && args+=(-v "$old_config:/old/config")
        volume_exists "$old_gh" && args+=(-v "$old_gh:/old/gh")
        args+=(-v "$repo/.devcontainer:/out")
        run_or_hint run "${args[@]}" alpine tar czf /out/.adopt.tgz --numeric-owner -C /old .
    fi

    # Unpacked by the container that already holds the volume. Copying only what
    # is missing is what keeps a credential from a fresh login; `cp -n` and
    # `tar --skip-old-files` are both absent from busybox, so it is spelled out.
    adevcontainer exec --name "$into" -- sh -c '
        set -e
        tmp="$(mktemp -d)"
        trap '"'"'rm -rf "$tmp"'"'"' EXIT
        tar xzf '"$stage_guest"' -C "$tmp"
        copy_missing() {
            [ -d "$1" ] || return 0
            ( cd "$1" && find . -type f ) | while read -r f; do
                [ -e "$2/$f" ] && continue
                mkdir -p "$2/$(dirname "$f")"
                cp -a "$1/$f" "$2/$f"
            done
        }
        mkdir -p /home/node/.state/claude /home/node/.state/gh
        copy_missing "$tmp/config" /home/node/.state/claude
        copy_missing "$tmp/gh" /home/node/.state/gh
        if [ -f "$tmp/hist/.bash_history" ] && [ ! -f /home/node/.state/history ]; then
            cp -a "$tmp/hist/.bash_history" /home/node/.state/history
        fi
        touch /home/node/.state/.adopted
    '

    rm -f "$stage_host"
    echo "state.sh: adopted. Anything already in the volume was left alone."
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

    # Staged inside the repository, not in `mktemp -d`. On macOS that is
    # /var/folders/…, and Apple `container` would not attach it: the bootstrap
    # fails with `VZErrorDomain Code=2 "The storage device attachment is
    # invalid."`, which names the storage layer and not the path, so it reads as
    # the volume-already-attached case and is not. `export` and `import` have
    # always bound a directory inside the repo, which is the shareable one; this
    # now does the same.
    #
    # 700 and removed on the way out, including on failure: it holds the Claude
    # credential for the length of one copy.
    staging="$repo/.state-adopt"
    rm -rf "$staging"
    mkdir -p "$staging"
    chmod 700 "$staging"
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

    # The volume attaches to one container at a time, so this project's own
    # container has to let go of it first. Stopping rather than asking you to:
    # adopt runs before you work in the container, so there is nothing in there
    # to interrupt, and `up.sh` starts it again. Quiet and ignored when it is not
    # running, which is the usual case.
    if "$rt" stop "$base" >/dev/null 2>&1; then
        echo "state.sh: stopped $base to free the volume; up.sh will start it again."
    fi

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
