#!/usr/bin/env bash
# Open this repository's devcontainer as a cmux *remote workspace* rather than a
# pane with a shell in it.
#
# cmux has no devcontainer support and does not need any: `cmux ssh` is a
# first-class workspace type with a Linux-side daemon, and a container running
# sshd is a remote host like any other.
#
# What it buys over `up.sh`, which is still the right tool for a quick shell:
#
#   - the `cmux` CLI works *inside* the container, so `cmux notify` and
#     `cmux workspace status set` from a Claude Code hook reach the app
#   - terminals survive cmux quitting, and reconnect on relaunch
#   - sidebar metadata, and dragging a file into the pane uploads over sftp
#
# Getting the container up and reachable is `ssh-target.sh`, shared with
# `orca-target.sh`. What is left here is the part that is about cmux.
#
#   ./.devcontainer/cmux-attach.sh             a workspace with a shell
#   ./.devcontainer/cmux-attach.sh claude      a workspace running Claude Code
#   REBUILD=1 ./.devcontainer/cmux-attach.sh   rebuild first
set -euo pipefail

# Called by `ssh-target.sh` after its "you are on Linux" guard and before it
# builds anything, which is the only correct place for it. Defining it here
# rather than checking inline above: checked above, this fires *first*, and on
# Linux -- inside a container, where cmux's CLI is absent until a `cmux ssh`
# session relays it -- the answer to "I am in the wrong shell" becomes "install
# cmux with brew", on Linux, where that goes nowhere. The guard has to win.
st_frontend_check() {
    command -v cmux >/dev/null 2>&1 && return 0
    echo "cmux-attach: cmux is not on PATH." >&2
    # cmux puts its CLI on the PATH of terminals *it* spawns, not on the system
    # one, so "not on PATH" usually means a plain Terminal.app window rather than
    # a missing install. Suggesting `brew install` first sends you to reinstall
    # something you already have.
    echo "cmux-attach:   run this from a cmux tab -- cmux only puts" >&2
    echo "cmux-attach:   its CLI on the PATH of terminals it starts." >&2
    echo "cmux-attach:   Not installed at all? brew install --cask cmux" >&2
    echo "cmux-attach: on Docker, use ./.devcontainer/up.sh instead." >&2
    exit 1
}

# Sourced without arguments so "$@" stays this script's -- it is the command to
# run in the workspace, and the helper has no business seeing it.
# shellcheck source=./ssh-target.sh
. "$(dirname "${BASH_SOURCE[0]}")/ssh-target.sh"

# sshd builds its own login environment, so what is on PATH there depends on
# image plumbing this script has no business assuming. npm's global bin -- where
# `claude` is -- is added here, single-quoted so $PATH expands in the container
# and not on the Mac. Belongs in the image too, and is, but a script that only
# works against a freshly built image is a script that fails at the worst time.
prelude='export PATH="$PATH:/usr/local/share/npm-global/bin"'

remote_command="$prelude && cd /workspace"
if [ "$#" -gt 0 ]; then
    remote_command="$prelude && cd /workspace && $*"
fi

# Pin how to get back here, from inside, where the relayed cmux CLI lives.
#
# cmux reconnects a saved ssh workspace on relaunch, and reconnecting is all it
# can do -- Apple `container` has no restart policy, so a container stopped while
# cmux was closed stays stopped and the workspace comes back as a remote daemon
# error. The only thing that recovers it is this script, which starts the
# container before attaching. A resume command is how cmux is told to run it.
#
# Set on every attach rather than once by hand, because "once by hand, per
# project, remembered" is the kind of step that is never done for the fourth
# project.
#
# Guarded rather than escaped: a quote in the path would need careful nesting
# through two shells, and a checked skip is worth more than clever quoting that
# is wrong once.
case "$repo" in
    *"'"*)
        echo "cmux-attach: path contains a quote; not pinning a resume command." >&2
        ;;
    *)
        # Double quotes around the path, single around the assignment. The
        # reverse -- which is what this said first -- ends the outer quoting at
        # the path and silently leaves it bare, which works until a path has a
        # space in it and then does not.
        resume="cd \"$repo\" && ./.devcontainer/cmux-attach.sh${*:+ $*}"
        # Manual restore is all this buys, so the message says so. It used to
        # say "approve it under Settings > Terminal > Resume Commands", which
        # sent the reader to a pane that stays at 0 commands forever: cmux
        # writes no approval record for a remote surface, and raises no prompt
        # for a CLI-set binding. A managed cmux ssh workspace is meant to come
        # back by remote PTY reattach instead, which is why the guards are
        # there. README.md, "What still does not work", cites the source.
        #
        # So this pin is not load-bearing. It is here so that a later session
        # can ask `cmux surface resume show` what a pane was for.
        #
        # Reported rather than silenced because while that was unknown there was
        # no way to tell a failed call from a missing cmux from one that never
        # ran. A line of output is cheaper than that question.
        pin='if command -v cmux >/dev/null 2>&1; then
    # --cwd explicitly. It defaults to $PWD, and $PWD in here is /workspace --
    # a path that does not exist on the Mac, where the restore actually runs, so
    # the binding would be stored pointing at nothing. No apostrophes in this
    # comment: the whole block is a single-quoted string and one ends it.
    if cmux surface resume set --cwd "$CMUX_ATTACH_CWD" --name "$CMUX_ATTACH_NAME" \
            --shell "$CMUX_ATTACH_RESUME" >/dev/null 2>&1; then
        echo "cmux-attach: resume command pinned for manual restore"
    else
        echo "cmux-attach: cmux surface resume set failed; nothing pinned" >&2
    fi
else
    echo "cmux-attach: no cmux CLI in the container, so nothing was pinned." >&2
    echo "cmux-attach: that means this is not a cmux ssh session -- up.sh does not provide one." >&2
fi'
        remote_command="$pin
$remote_command"
        remote_command="CMUX_ATTACH_RESUME='$resume'; CMUX_ATTACH_CWD='$repo'; CMUX_ATTACH_NAME='$container'; $remote_command"
        ;;
esac

# Host key checking off: start-sshd.sh generates a host key per container, so
# known_hosts could only ever reject a rebuild of the same workspace. What bounds
# this is the address being the runtime's own and the key being one this made.
#
# --no-forward-agent for the reason the container holds its own gh credential:
# it is a sandbox, and the host's ssh agent is not part of what it gets.
exec cmux ssh "node@$target" \
    --port 2222 \
    --identity "$key" \
    --name "$container" \
    --no-forward-agent \
    --ssh-option IdentitiesOnly=yes \
    --ssh-option UserKnownHostsFile=/dev/null \
    --ssh-option StrictHostKeyChecking=no \
    --command "$remote_command"
