#!/usr/bin/env bash
# A daily snapshot of this container's state into somewhere durable, run by
# launchd while you are already using the Mac.
#
#   ./.devcontainer/state-backup.sh --install ~/OneDrive/devcontainer-state
#   ./.devcontainer/state-backup.sh --run          once, by hand
#   ./.devcontainer/state-backup.sh --status
#   ./.devcontainer/state-backup.sh --uninstall
#
# `StartInterval`, deliberately, and not `StartCalendarInterval`: a calendar
# entry asks to happen at a wall-clock time and can bring the machine up to
# honour it. An interval only says "no more than a day between runs", so a
# sleeping Mac stays asleep and the run happens shortly after it next wakes.
# `LimitLoadToSessionType Aqua` adds the other half: it loads for a logged-in
# GUI session, so this only exists while you are actually using the machine.
#
# The volume is not synced and must not be -- it is a live filesystem, and a
# sync engine copying one while it is mounted produces torn writes. A tarball is
# a consistent snapshot and an ordinary file, which is the thing cloud storage
# is good at. See README.md.
set -euo pipefail

if [ "$(uname -s)" = Linux ]; then
    printf '%s\n' \
        "${0##*/}: this runs on the Mac, not inside the container." \
        "${0##*/}: \`exit\` back to the host first, or use another cmux tab." >&2
    exit 1
fi

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
base="$(basename "$repo")"
label="com.orhunkok.devcontainer-state.$base"
plist="$HOME/Library/LaunchAgents/$label.plist"
log="$HOME/Library/Logs/$label.log"

usage() { sed -n '2,9p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

# --- the run itself -------------------------------------------------------

do_run() {
    dest="${STATE_BACKUP_DEST:-${1:-}}"
    if [ -z "$dest" ]; then
        echo "state-backup: no destination. Pass one, or set STATE_BACKUP_DEST." >&2
        exit 1
    fi

    # Quietly do nothing rather than fail: an agent that reports an error every
    # day because the runtime happens to be down is an agent you stop reading.
    if ! command -v container >/dev/null 2>&1 && ! command -v docker >/dev/null 2>&1; then
        echo "$(date '+%F %T') no container runtime; skipping."
        exit 0
    fi
    if ! container system status >/dev/null 2>&1 && ! docker info >/dev/null 2>&1; then
        echo "$(date '+%F %T') no runtime running; skipping."
        exit 0
    fi

    mkdir -p "$dest"
    archive="$dest/$base.tar.gz"

    # Two deep. A single overwritten file is a backup that a bad export destroys,
    # and a bad export is exactly what you would not notice for a month.
    [ -f "$archive" ] && mv -f "$archive" "$dest/$base.prev.tar.gz"

    # No --with-credentials: this lands in cloud storage, and a GitHub token
    # does not belong there. `gh auth login` is one command on the far side.
    if "$repo/.devcontainer/state.sh" export "$archive"; then
        echo "$(date '+%F %T') wrote $archive ($(du -h "$archive" | cut -f1))"
    else
        echo "$(date '+%F %T') export failed; keeping the previous snapshot." >&2
        [ -f "$dest/$base.prev.tar.gz" ] && mv -f "$dest/$base.prev.tar.gz" "$archive"
        exit 1
    fi
}

# --- install / remove ----------------------------------------------------

do_install() {
    dest="${1:-}"
    [ -n "$dest" ] || { echo "state-backup: --install needs a destination directory." >&2; exit 1; }
    dest="$(cd "$(dirname "$dest")" && pwd)/$(basename "$dest")"

    mkdir -p "$(dirname "$plist")" "$(dirname "$log")" "$dest"

    # launchd hands a process a minimal PATH, so `container` and anything from
    # Homebrew have to be named. This is the usual reason a working script fails
    # only under launchd.
    cat > "$plist" <<PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>Label</key><string>$label</string>
    <key>ProgramArguments</key>
    <array>
        <string>/bin/bash</string>
        <string>$repo/.devcontainer/state-backup.sh</string>
        <string>--run</string>
    </array>
    <key>StartInterval</key><integer>86400</integer>
    <key>RunAtLoad</key><false/>
    <key>LimitLoadToSessionType</key><string>Aqua</string>
    <key>ProcessType</key><string>Background</string>
    <key>LowPriorityIO</key><true/>
    <key>Nice</key><integer>5</integer>
    <key>WorkingDirectory</key><string>$repo</string>
    <key>StandardOutPath</key><string>$log</string>
    <key>StandardErrorPath</key><string>$log</string>
    <key>EnvironmentVariables</key>
    <dict>
        <key>PATH</key><string>/opt/homebrew/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin</string>
        <key>STATE_BACKUP_DEST</key><string>$dest</string>
    </dict>
</dict>
</plist>
PLIST

    plutil -lint "$plist" >/dev/null

    # bootout first so --install is idempotent; both are allowed to fail, since
    # the first run has nothing loaded to remove.
    launchctl bootout "gui/$(id -u)/$label" 2>/dev/null || true
    launchctl bootstrap "gui/$(id -u)" "$plist" 2>/dev/null \
        || launchctl load "$plist"

    echo "state-backup: installed $label"
    echo "state-backup:   snapshots -> $dest/$base.tar.gz (plus .prev)"
    echo "state-backup:   log       -> $log"
    echo "state-backup:   at most a day apart, only while logged in, never wakes the Mac"
    echo "state-backup: run one now with: $0 --run $dest"
}

do_uninstall() {
    launchctl bootout "gui/$(id -u)/$label" 2>/dev/null \
        || launchctl unload "$plist" 2>/dev/null \
        || true
    rm -f "$plist"
    echo "state-backup: removed $label (snapshots and log left alone)"
}

do_status() {
    if launchctl print "gui/$(id -u)/$label" >/dev/null 2>&1; then
        echo "state-backup: $label is loaded"
    else
        echo "state-backup: $label is not loaded"
    fi
    [ -f "$plist" ] && echo "state-backup: plist $plist" || echo "state-backup: no plist"
    [ -f "$log" ] && tail -3 "$log" || echo "state-backup: no log yet"
}

case "${1:-}" in
    --run) shift; do_run "$@" ;;
    --install) shift; do_install "$@" ;;
    --uninstall) do_uninstall ;;
    --status) do_status ;;
    -h|--help|"") usage ;;
    *) echo "state-backup: unknown option $1" >&2; usage >&2; exit 1 ;;
esac
