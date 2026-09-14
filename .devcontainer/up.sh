#!/usr/bin/env bash
# Bring this repository's devcontainer up and open a shell in it -- the job the
# VS Code Dev Containers extension used to do, minus the editor.
#
# This runs on the host, not inside the container. Call it by hand, or from
# whatever your terminal app uses to pin a command to a pane; see README.md.
#
#   ./.devcontainer/up.sh             a zsh inside the container
#   ./.devcontainer/up.sh claude      straight into Claude Code
#   REBUILD=1 ./.devcontainer/up.sh   discard the existing container first
#
# Apple `container`, through `adevcontainer`. The config is image-based, so
# `build.sh` makes the image and has to have been run at least once.
set -euo pipefail

# These drive macOS-side tooling -- `container`, `adevcontainer` -- so running one
# *inside* the container is a mistake worth naming. Left uncaught the symptom is a
# tool "not on PATH" plus an invitation to `brew install` it, on Linux, which
# sends you somewhere with no exit.
if [ "$(uname -s)" = Linux ]; then
    printf '%s\n' \
        "${0##*/}: this runs on the Mac, not inside the container." \
        "${0##*/}: \`exit\` back to the host first, or open a terminal there." >&2
    exit 1
fi

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if ! command -v adevcontainer >/dev/null 2>&1; then
    echo "up.sh: adevcontainer is not on PATH." >&2
    echo "up.sh:   brew install wcgomes/tap/adevcontainer" >&2
    exit 1
fi

# adevcontainer discovers .devcontainer/devcontainer.json from the working
# directory and takes no --workspace-folder.
cd "$repo"

# One log for both branches below, and it has to outlive them: the container name
# is read out of it. `rebuild` previously wrote no log at all and `up` deleted its
# own before anything could read it, so the name always fell back to the
# directory -- which is right only while the folder and `name` in
# devcontainer.json agree.
log="$(mktemp)"
trap 'rm -f "$log"' EXIT

if [ -n "${REBUILD:-}" ]; then
    # `rebuild --name` selects an *existing* container, so the name has to be the
    # one that exists -- not the one devcontainer.json would create. Those differ
    # whenever `name` has been edited or the folder renamed, and the failure is
    # `No managed container named ...` while the container is sitting right there
    # under another name. `list` knows which container belongs to this folder;
    # ask it. The leading slash in the match keeps `notpinpoint` from answering
    # for `pinpoint`.
    existing="$(adevcontainer list 2>/dev/null \
        | awk -v b="/$(basename "$repo")" \
            'NR>1 { n=length(b); if (substr($NF, length($NF)-n+1) == b) { print $1; exit } }' \
    || true)"

    if [ -n "$existing" ]; then
        adevcontainer rebuild --name "$existing" 2>&1 | tee "$log"
    else
        echo "up.sh: no container for this folder yet; creating one." >&2
        adevcontainer up 2>&1 | tee "$log"
    fi
else
    # `up` fails closed when devcontainer.json has changed since the container
    # was created, which is correct of it and a dead end here: its hint names
    # `adevcontainer rebuild` rather than the way in from this script. Translate
    # rather than make you map it. Not automatic, because a rebuild replaces the
    # container and would take an agent running in another pane with it --
    # volumes are preserved, work in progress is not.
    #
    # Streamed through `tee` rather than captured into a variable. This read the
    # output into `$out` and printed it afterwards, which meant a first start --
    # image pull, VM boot, postCreate -- showed nothing at all for minutes and
    # was indistinguishable from a hang. `pipefail` is set, so the pipeline still
    # fails when adevcontainer does, and the log is kept only to grep for the one
    # error worth translating.
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

# --name for the same reason ssh-target.sh uses it: with two managed containers
# running, `exec` without it opens an interactive picker and acts on whichever row
# is highlighted. Taken from the tool's own `containerId:`.
container="$(sed -n 's/.*containerId:[[:space:]]*\([A-Za-z0-9_.-][A-Za-z0-9_.-]*\).*/\1/p' "$log" | tail -1)"
rm -f "$log"
trap - EXIT   # this script ends in exec, which would never reach the trap
[ -n "$container" ] || container="$(basename "$repo")"

# Once per volume, bring this project's earlier Claude state in if it has any.
# The marker check is one cheap exec against a container that is already running;
# the search behind it only happens the first time. Reading a *Docker* volume is
# part of that search, which is why state.sh still knows about Docker when
# nothing else here does -- see its header.
if ! adevcontainer exec --name "$container" -- test -e /home/node/.state/.adopted 2>/dev/null; then
    "$repo/.devcontainer/state.sh" adopt || \
        echo "up.sh: could not adopt earlier state; continuing." >&2
fi

# Nothing is passed with an env override here, and one measured thing is worth
# keeping in case it ever has to be. `docker exec -t` defaults TERM to plain
# `xterm` -- eight colours by terminfo -- and Claude Code's TUI drops to sixteen
# without saying so; the Docker path this script used to carry pinned
# `TERM=xterm-256color` against exactly that. `adevcontainer exec` has not been
# measured for the same fault, so nothing is pinned for it rather than pinning
# blind: a TERM that does not resolve breaks more than sixteen colours does. The
# ssh path is immune either way, because ssh carries TERM itself and
# sshd-remote.conf lists COLORTERM in AcceptEnv.
#
# `zsh`, not `zsh -l`: it matches the profile VS Code's terminal used, and it is
# the interactive shell that sources ~/.zshrc, which is where gh-auth.sh lives.
exec adevcontainer exec -it --name "$container" -- "${@:-zsh}"
