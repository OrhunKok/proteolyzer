#!/usr/bin/env bash
# Make the cmux sidebar say what the agent in this container is doing.
#
# Run inside the container, once per ~/.claude volume:
#
#   /workspace/.devcontainer/install-cmux-hooks.sh
#
# Why this is needed at all: cmux only draws its agent status pill for a process
# it recognises as an agent, and it recognises the `claude` it wrapped on the
# host -- not one running inside a container, which it sees as a generic
# process. AgentBox hit the same wall and worked around it by driving the
# workspace's own presentation instead. `cmux ssh` puts the CLI on PATH in here
# with CMUX_SOCKET_PATH pinned to this workspace's relay, so the hooks below
# reach the right workspace with no host-side daemon in the middle.
#
# `cmux notify`, and not `cmux workspace status set`, which is what this wrote
# until 2026-09-11 and never worked. Inside the container `cmux` is a wrapper
# that execs `cmuxd-remote`, a different binary from the host CLI: its only
# `workspace` subcommand is `group`, so the status call was an argument parse
# error -- exit 2, before it ever reached the relay. `|| true` swallowed it, so
# there was no error and the pill simply never moved. The daemon's advertised
# capabilities (`cmux capabilities`) carry no status lane at all, and
# `needs-attention` appears nowhere in the binary, so this is not a rename
# waiting to be found.
#
# `notify` is the lane that exists: it takes no positional arguments, reads the
# surface and workspace from the cmux env vars, and badges the workspace in the
# sidebar. Verified end to end on 2026-09-11 -- `cmux notify` then
# `cmux jump-to-unread` returns the notification it made, with the body below.
#
# Nothing happens on the `up.sh` path: there is no `cmux` on PATH there, and the
# guard in each hook makes it a silent no-op rather than an error every turn.
set -euo pipefail

settings="${CLAUDE_CONFIG_DIR:-$HOME/.claude}/settings.json"

if ! command -v jq >/dev/null 2>&1; then
    echo "install-cmux-hooks: jq is required and is in this image; are you on the host?" >&2
    exit 1
fi

# `|| true` so a hook can never fail a turn, and the `command -v` guard so this
# is inert outside a cmux ssh workspace. `timeout` because a hook that blocks
# blocks the turn: a dead relay refuses instantly, but a wedged one would not,
# and CMUX_SOCKET_PATH has been seen pointing at a port nothing was listening on
# after a reconnect.
#
# Both hooks notify, with different bodies, because both are the same event from
# the sidebar's point of view: the agent has stopped and it is your turn. Stop
# fires when a response completes, which for a container you left running is the
# signal worth having. If that is too chatty while you are watching the pane,
# `cmux mark-notification-read` is the one-word swap for the Stop line -- it
# clears rather than raises, and takes no arguments either.
notify='command -v cmux >/dev/null 2>&1 && timeout 5 cmux notify --title Claude --body "needs your input" || true'
stop='command -v cmux >/dev/null 2>&1 && timeout 5 cmux notify --title Claude --body "finished its turn" || true'

mkdir -p "$(dirname "$settings")"
[ -f "$settings" ] || echo '{}' > "$settings"

# Idempotent: drop any entry this script wrote before -- matched on the command
# text -- then append the current pair. Hooks from anywhere else are preserved.
# The match covers the dead `workspace status` spelling so an old install is
# cleaned out rather than left beside the new one, and `mark-notification-read`
# so the Stop swap above stays re-runnable too.
tmp="$(mktemp)"
jq --arg notify "$notify" --arg stop "$stop" '
  def mine: tostring
    | (contains("cmux workspace status")
       or contains("cmux notify")
       or contains("cmux mark-notification-read"));
  (.hooks // {}) as $h
  | ($h.Notification // [] | map(select(mine | not))) as $keptNotification
  | ($h.Stop         // [] | map(select(mine | not))) as $keptStop
  | .hooks = ($h
      | .Notification = ($keptNotification + [{hooks: [{type: "command", command: $notify}]}])
      | .Stop         = ($keptStop         + [{hooks: [{type: "command", command: $stop}]}]))
' "$settings" > "$tmp"
mv "$tmp" "$settings"

echo "install-cmux-hooks: wrote Notification and Stop hooks to $settings"
echo "install-cmux-hooks: the workspace badges in the sidebar when Claude hands back."
