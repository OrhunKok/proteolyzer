#!/bin/bash
# Start the sshd that `cmux ssh` attaches to. Idempotent, run from
# postStartCommand, and reachable through sudo for the same reason
# init-firewall.sh is: it needs root and nothing else here does.
#
# Starting it unconditionally is safe. The port is published on loopback only,
# password auth is off, and no authorized_keys exists until cmux-attach.sh
# writes one -- so until you ask for it, this listens and refuses everything.
set -euo pipefail

CONFIG=/etc/ssh/sshd_cmux_config
PIDFILE=/run/sshd-cmux.pid

# Generated at first start rather than baked into the image, so two containers
# from the same image do not share a host key.
if [ ! -f /etc/ssh/ssh_host_ed25519_key ]; then
    ssh-keygen -A
fi

# Debian's sshd wants its privilege separation directory to exist.
mkdir -p /run/sshd

if [ -f "$PIDFILE" ] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
    echo "start-sshd: already running on port 2222 (pid $(cat "$PIDFILE"))"
    exit 0
fi

# `sshd -t` does not check Subsystem paths, and that one is distribution-
# specific: /usr/lib/openssh on Debian, /usr/lib/ssh on Alpine. Only the sftp
# upload path cares, so say so rather than refusing to start.
if [ ! -x /usr/lib/openssh/sftp-server ]; then
    echo "start-sshd: warning: /usr/lib/openssh/sftp-server missing;" >&2
    echo "start-sshd: dragging a file into a cmux pane will not upload." >&2
fi

# Rejects its own config before backgrounding, so a mistake in the file shows up
# here rather than as a connection refused ten minutes later.
/usr/sbin/sshd -t -f "$CONFIG"
/usr/sbin/sshd -f "$CONFIG"
echo "start-sshd: listening on port 2222"
