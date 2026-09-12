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
#
# Uses `adevcontainer` (Apple `container`) when it is installed and the
# devcontainer CLI (Docker) otherwise. The config is image-based, so both read
# it; `build.sh` makes the image and has to have been run at least once.
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

if command -v adevcontainer >/dev/null 2>&1; then
    cli=adevcontainer
elif command -v devcontainer >/dev/null 2>&1; then
    cli=devcontainer
else
    echo "up.sh: no devcontainer CLI found." >&2
    echo "up.sh:   brew install wcgomes/tap/adevcontainer   (Apple container)" >&2
    echo "up.sh:   npm install -g @devcontainers/cli        (Docker)" >&2
    exit 1
fi

if [ "$cli" = adevcontainer ]; then
    # adevcontainer discovers .devcontainer/devcontainer.json from the working
    # directory and takes no --workspace-folder.
    cd "$repo"

    # One log for both branches, and it has to outlive them: the container name
    # is read out of it below. `rebuild` previously wrote no log at all and `up`
    # deleted its own before anything could read it, so the name always fell
    # back to the directory -- which is right only while the folder and `name`
    # in devcontainer.json agree.
    log="$(mktemp)"
    trap 'rm -f "$log"' EXIT

    if [ -n "${REBUILD:-}" ]; then
        adevcontainer rebuild 2>&1 | tee "$log"
    else
        # `up` fails closed when devcontainer.json has changed since the
        # container was created, which is correct of it and a dead end here: its
        # hint names `adevcontainer rebuild` rather than the way in from this
        # script. Translate rather than make you map it. Not automatic, because
        # a rebuild replaces the container and would take an agent running in
        # another pane with it -- volumes are preserved, work in progress is not.
        #
        # Streamed through `tee` rather than captured into a variable. This read
        # the output into `$out` and printed it afterwards, which meant a first
        # start -- image pull, VM boot, postCreate -- showed nothing at all for
        # minutes and was indistinguishable from a hang. The Docker branch below
        # already says exactly that about itself; this branch was the one doing
        # it. `pipefail` is set, so the pipeline still fails when adevcontainer
        # does, and the log is kept only to grep for the one error worth
        # translating.
        if ! adevcontainer up 2>&1 | tee "$log"; then
            if grep -q config_hash "$log"; then
                echo >&2
                echo "up.sh: devcontainer.json changed since this container was made." >&2
                echo "up.sh: rebuild it -- volumes are kept, the container is replaced:" >&2
                echo "up.sh:   REBUILD=1 $0${*:+ $*}" >&2
            fi
            exit 1
        fi
    fi
    # --name for the same reason cmux-attach.sh uses it: with two managed
    # containers running, `exec` without it opens an interactive picker and acts
    # on whichever row is highlighted. Taken from the tool's own `containerId:`.
    container="$(sed -n 's/.*containerId:[[:space:]]*\([A-Za-z0-9_.-][A-Za-z0-9_.-]*\).*/\1/p' "$log" | tail -1)"
    rm -f "$log"
    trap - EXIT   # this branch ends in exec, which would never reach the trap
    [ -n "$container" ] || container="$(basename "$repo")"
    exec adevcontainer exec -it --name "$container" -- "${@:-zsh}"
fi

if ! docker info >/dev/null 2>&1; then
    echo "up.sh: Docker is not running." >&2
    exit 1
fi

# Idempotent: reuses a container that is already up. Output is left visible
# because a first start takes a while and silence there reads as a hang.
# shellcheck disable=SC2086
devcontainer up --workspace-folder "$repo" ${REBUILD:+--remove-existing-container}

# cmux sets the CMUX_ ones in its own terminals, and they are how anything inside
# the container says which pane it is talking about -- a notification from a hook
# is otherwise untargeted. COLORTERM rides along for the reason below. Only
# forwarded when set, so running this outside cmux does not plant empty
# variables. The idea is lifted from zackey-heuristics/cmux-devcontainer-bridge,
# which needs the same thing.
remote_env=()
for var in CMUX_WORKSPACE_ID CMUX_SURFACE_ID CMUX_TAB_ID COLORTERM; do
    value="${!var:-}"
    if [ -n "$value" ]; then
        remote_env+=(--remote-env "$var=$value")
    fi
done

# `docker exec -t` defaults TERM to plain `xterm`, which terminfo says is eight
# colours, and Claude Code's TUI drops to sixteen and says nothing about it. The
# ssh path never had this -- ssh carries TERM itself and sshd-cmux.conf lists
# COLORTERM in AcceptEnv -- so it is the exec path alone that arrives washed out.
#
# Pinned, not forwarded from the host: the container's terminfo is Debian's, and
# a host running Ghostty or kitty exports an entry (xterm-ghostty, xterm-kitty)
# that is not in it. A TERM that does not resolve breaks considerably more than
# sixteen colours does, and xterm-256color is both present here and what those
# terminals fall back to anyway.
remote_env+=(--remote-env "TERM=xterm-256color")

# `zsh`, not `zsh -l`: it matches the profile VS Code's terminal used, and it is
# the interactive shell that sources ~/.zshrc, which is where gh-auth.sh lives.
# The ${x[@]+"${x[@]}"} form is for macOS's bash 3.2, where expanding an empty
# array under `set -u` is an error.
exec devcontainer exec --workspace-folder "$repo" \
    ${remote_env[@]+"${remote_env[@]}"} \
    "${@:-zsh}"
