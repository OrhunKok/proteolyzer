#!/usr/bin/env bash
# Make this repository's devcontainer ready for Orca, and print what to type into
# Settings -> SSH -> Add Target.
#
# Orca (stablyai/orca) is the frontend this container is served for. It has no
# devcontainer support and does not need any: an SSH target with a custom port
# and identity file is a first-class thing in its UI, and a container running
# sshd is a remote host like any other.
#
# What that buys, and why this file exists:
#
#   - a file tree, a fuzzy finder, an editor and a diff view over the same ssh
#     connection, operating on the container's filesystem
#   - remote PTYs leased by a relay *on the container*, so they survive Orca
#     quitting and come back attached with their scrollback. See README.md.
#
# There is no CLI to hand a target to -- targets are added in the app. So this
# does the parts a script can do (start the container, put the key in, start
# sshd, prove the login, prove Orca's relay will build) and then prints the four
# fields, rather than pretending to automate a GUI form.
#
#   ./.devcontainer/orca-target.sh             prepare and print
#   REBUILD=1 ./.devcontainer/orca-target.sh   rebuild first
set -euo pipefail

# shellcheck source=./ssh-target.sh
. "$(dirname "${BASH_SOURCE[0]}")/ssh-target.sh"

ssh_opts=(
    -o BatchMode=yes
    -o IdentitiesOnly=yes
    -o UserKnownHostsFile=/dev/null
    -o StrictHostKeyChecking=no
    -o ConnectTimeout=10
    -i "$key" -p 2222
)

# The check that decides whether this is worth doing at all.
#
# On first connect Orca installs a relay on the remote, and its remote terminals
# need a native `node-pty`. Linux has no prebuild, so it compiles on the host and
# wants make, a C++ compiler and python3. Orca's own docs are explicit that
# without them connect still succeeds for "files, git, and the editor" and
# terminals simply do not work -- which, for running an agent, is the whole
# point. A green file tree and a dead terminal is a confusing way to find out.
#
# Deliberately over `ssh <host> '<command>'` and not an interactive shell. That
# is a non-login, non-interactive shell, which is what Orca's installer gets, and
# it is the environment where things go missing: container `ENV` does not reach
# an sshd session, and this image's /etc/profile.d and /etc/zsh/zprofile wiring
# only runs for *login* shells. `claude not found` was this same gap. What saves
# it is sshd's built-in default PATH containing /usr/local/bin, which is why the
# image symlinks what has to be found into there.
echo
echo "orca-target: checking what Orca's relay needs, the way its installer sees it"
missing=""
for tool in make g++ python3; do
    if found="$(ssh "${ssh_opts[@]}" "node@$target" "command -v $tool" 2>/dev/null)" \
            && [ -n "$found" ]; then
        printf '  %-8s %s\n' "$tool" "$found"
    else
        printf '  %-8s MISSING\n' "$tool"
        missing="$missing $tool"
    fi
done

if [ -n "$missing" ]; then
    echo
    echo "orca-target: remote terminals will NOT work:$missing not found on a" >&2
    echo "orca-target: non-interactive ssh PATH. Files, git and the editor still" >&2
    echo "orca-target: will. The fix is in the image, not here -- symlink the" >&2
    echo "orca-target: missing tool into /usr/local/bin, which is on sshd's" >&2
    echo "orca-target: default PATH, then REBUILD=1 and run this again." >&2
    echo
else
    echo "  all three found -- node-pty will build, so remote terminals will work"
fi

# A relative base path, and the reason is not style.
#
# For an SSH repo Orca ignores an absolute global workspace dir and falls back to
# a sibling of the repo: `join(repoPath, '..', repoName + '-' + name)` in
# src/main/ipc/worktree-logic.ts. The repo here is /workspace, so a sibling is
# /<name> at the *container root* -- container-local disk, invisible from the
# Mac, outside the state volume, and gone on the next rebuild. A relative value
# resolves inside the repo instead, which is the Mac's disk over virtiofs and
# survives everything.
#
# .claude/worktrees is what this repository already uses for its own worktrees,
# so Orca's land beside them rather than inventing a second convention.
worktrees=".claude/worktrees"

cat <<EOF

orca-target: Settings -> SSH -> Add Target

  Host            $target
  User            node
  Port            2222
  Identity file   $key

Then "Test", then "Save". Add the repo with this target as its location and the
remote path below, or open the folder directly from the file picker.

  Remote path     /workspace

Then one setting, in the repo's worktree base path:

  Worktree base   $worktrees

It has to be *relative*. Orca ignores an absolute base path for an ssh repo and
puts worktrees in a sibling of the repo instead -- and the repo here is
/workspace, so a sibling is the container's root filesystem: invisible from the
Mac, outside the state volume, and destroyed by the next rebuild. A relative
path lands inside /workspace, which is the Mac's disk.

EOF

if [ "$target" = "$ip" ]; then
    cat >&2 <<EOF
orca-target: this is the address, not a name, because this Mac does not resolve
orca-target:   $host
orca-target: An Orca target is a saved host entry and the address changes on
orca-target: every rebuild, so the target will go stale and need editing. Set the
orca-target: DNS domain up and it stops mattering:
orca-target:   container system dns create ${CONTAINER_DNS_DOMAIN:-adevcontainers.local}

EOF
fi

echo "orca-target: the sshd this reaches is not Orca-specific -- anything that"
echo "orca-target: takes a host, a port and an identity file can read the same one."
