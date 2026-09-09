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
# workspace's own presentation instead. cmux exposes a documented lane for
# exactly that, `cmux workspace status set`, and `cmux ssh` puts the CLI on PATH
# in here with CMUX_SOCKET_PATH pinned to this workspace's relay. So the hooks
# below reach the right workspace with no host-side daemon in the middle.
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
# is inert outside a cmux ssh workspace.
notify='command -v cmux >/dev/null 2>&1 && cmux workspace status set needs-attention || true'
stop='command -v cmux >/dev/null 2>&1 && cmux workspace status set auto || true'

mkdir -p "$(dirname "$settings")"
[ -f "$settings" ] || echo '{}' > "$settings"

# Idempotent: drop any entry this script wrote before -- matched on the command
# text -- then append the current pair. Hooks from anywhere else are preserved.
tmp="$(mktemp)"
jq --arg notify "$notify" --arg stop "$stop" '
  (.hooks // {}) as $h
  | ($h.Notification // [] | map(select(tostring | contains("cmux workspace status") | not))) as $keptNotification
  | ($h.Stop         // [] | map(select(tostring | contains("cmux workspace status") | not))) as $keptStop
  | .hooks = ($h
      | .Notification = ($keptNotification + [{hooks: [{type: "command", command: $notify}]}])
      | .Stop         = ($keptStop         + [{hooks: [{type: "command", command: $stop}]}]))
' "$settings" > "$tmp"
mv "$tmp" "$settings"

echo "install-cmux-hooks: wrote Notification and Stop hooks to $settings"
echo "install-cmux-hooks: the sidebar row turns amber when Claude wants you."
