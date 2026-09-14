FROM node:20

ARG TZ
ENV TZ="$TZ"

ARG CLAUDE_CODE_VERSION=latest

# Install basic development tools and iptables/ipset
RUN apt-get update && apt-get install -y --no-install-recommends \
  less \
  git \
  procps \
  sudo \
  fzf \
  zsh \
  man-db \
  unzip \
  gnupg2 \
  gh \
  iptables \
  ipset \
  iproute2 \
  dnsutils \
  aggregate \
  jq \
  nano \
  vim \
  openssh-server \
  && apt-get clean && rm -rf /var/lib/apt/lists/*

# Ensure default node user has access to /usr/local/share
RUN mkdir -p /usr/local/share/npm-global && \
  chown -R node:node /usr/local/share

ARG USERNAME=node

# Persist bash history.
RUN SNIPPET="export PROMPT_COMMAND='history -a' && export HISTFILE=/home/node/.state/history" \
  && mkdir -p /home/node/.state \
  && touch /home/node/.state/history \
  && chown -R $USERNAME /home/node/.state

# Set `DEVCONTAINER` environment variable to help with orientation
ENV DEVCONTAINER=true

# Create the workspace and the one state directory, node-owned.
#
# Everything that must survive a rebuild lives under /home/node/.state, which is
# the single named volume: `claude/` for Claude Code's config, projects and
# memory, `gh/` for the GitHub login, `git/` for the git identity, and `history`
# for the shell. There were three volumes and one is easier to keep track of;
# CLAUDE_CONFIG_DIR, GH_CONFIG_DIR and GIT_CONFIG_GLOBAL do the redirecting, and
# HISTFILE above.
#
# Every subdirectory is created here, node-owned, and that is load-bearing: a
# named volume is initialised from the image path it covers, ownership included,
# so a path the image does not have comes up root-owned and unwritable by node --
# with no warning, and sudo here is restricted to named scripts. decoder hit
# exactly that on a shared memory volume: five files failed to restore and the
# script died on a bare Permission denied. Apple `container` does not seed
# volumes at all, which fix-volume-perms.sh covers.
#
# git/config is seeded with an include rather than left empty, and that include
# is the whole reason this works. GIT_CONFIG_GLOBAL moves the global file onto
# the volume; ~/.gitconfig is where adevcontainer writes its credential helper at
# every start. Relocating the global file without pulling that back in would stop
# git finding the helper, and `git push` from inside the container would fail
# instead. A missing include path is ignored by git, so it costs nothing on a
# runtime that never writes one.
RUN mkdir -p /workspace \
  /home/node/.state/claude/projects/-workspace/memory \
  /home/node/.state/gh \
  /home/node/.state/git && \
  printf '%s\n' '[include]' '  path = /home/node/.gitconfig' \
    > /home/node/.state/git/config && \
  chown -R node:node /workspace /home/node/.state

WORKDIR /workspace

ARG GIT_DELTA_VERSION=0.18.2
RUN ARCH=$(dpkg --print-architecture) && \
  wget "https://github.com/dandavison/delta/releases/download/${GIT_DELTA_VERSION}/git-delta_${GIT_DELTA_VERSION}_${ARCH}.deb" && \
  sudo dpkg -i "git-delta_${GIT_DELTA_VERSION}_${ARCH}.deb" && \
  rm "git-delta_${GIT_DELTA_VERSION}_${ARCH}.deb"

# Set up non-root user
USER node

# Install global packages
ENV NPM_CONFIG_PREFIX=/usr/local/share/npm-global
ENV PATH=$PATH:/usr/local/share/npm-global/bin

# Set the default shell to zsh rather than sh
ENV SHELL=/bin/zsh

# Set the default editor and visual
ENV EDITOR=nano
ENV VISUAL=nano

# Default powerline10k theme
ARG ZSH_IN_DOCKER_VERSION=1.2.0
RUN sh -c "$(wget -O- https://github.com/deluan/zsh-in-docker/releases/download/v${ZSH_IN_DOCKER_VERSION}/zsh-in-docker.sh)" -- \
  -p git \
  -p fzf \
  -a "source /usr/share/doc/fzf/examples/key-bindings.zsh" \
  -a "source /usr/share/doc/fzf/examples/completion.zsh" \
  -a "export PROMPT_COMMAND='history -a' && export HISTFILE=/home/node/.state/history" \
  -x

# Install Claude
RUN npm install -g @anthropic-ai/claude-code@${CLAUDE_CODE_VERSION}

# The settings that used to live in devcontainer.json's containerEnv. They are
# here instead because a container's environment is given to processes the
# *runtime* starts -- an entrypoint, `docker exec`, `devcontainer exec` -- and
# sshd is not one of those. It builds a login environment from /etc/profile and
# friends, so over ssh none of these were set and `claude` was not even on PATH.
# One place, applied to both.
ENV CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC=true
ENV DISABLE_TELEMETRY=true
ENV NODE_OPTIONS=--max-old-space-size=4096
ENV CLAUDE_CONFIG_DIR=/home/node/.state/claude
ENV GH_CONFIG_DIR=/home/node/.state/gh
# git's global config, moved onto the volume for the same reason as those two:
# /home/node is not mounted, so an identity set with `git config --global` is
# gone on the next rebuild, and `git commit` then fails with "Author identity
# unknown" in a container that worked yesterday. The path is set here; the value
# is not, because publish.sh pushes this image to GHCR and an identity belongs to
# a person rather than to an image.
ENV GIT_CONFIG_GLOBAL=/home/node/.state/git/config
ENV POWERLEVEL9K_DISABLE_GITSTATUS=true

# The same values again for login shells, which read files rather than the
# container config.
#
# /etc/profile sources /etc/profile.d/*.sh, which covers bash. It does not cover
# zsh: this image's /etc/zsh/zprofile is the packaged one and contains nothing
# but comments, so a zsh login shell reads neither. sshd hands out zsh here
# (usermod --shell, below), which is precisely the shell that needed this, so
# zprofile gets the line Debian would otherwise have put there. `emulate sh`
# keeps /etc/profile's sh syntax from meeting zsh's parser.
USER root
RUN grep -q '/etc/profile' /etc/zsh/zprofile || \
  printf "%s\n" "emulate sh -c 'source /etc/profile'" >> /etc/zsh/zprofile
RUN printf '%s\n' \
  '# Written by .devcontainer/Dockerfile. A login shell -- notably one sshd' \
  '# starts -- does not inherit the container ENV, so it is restated here.' \
  'export PATH="$PATH:/usr/local/share/npm-global/bin"' \
  'export CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC=true' \
  'export DISABLE_TELEMETRY=true' \
  'export NODE_OPTIONS=--max-old-space-size=4096' \
  'export CLAUDE_CONFIG_DIR=/home/node/.state/claude' \
  'export GH_CONFIG_DIR=/home/node/.state/gh' \
  'export GIT_CONFIG_GLOBAL=/home/node/.state/git/config' \
  'export POWERLEVEL9K_DISABLE_GITSTATUS=true' \
  > /etc/profile.d/10-devcontainer.sh && \
  chmod 0644 /etc/profile.d/10-devcontainer.sh

# Belt and braces for a shell that reads no profile at all -- `bash -lc` resets
# PATH to the Debian default, and some agents shell out that way. /usr/local/bin
# is in every default PATH there is.
RUN ln -sf /usr/local/share/npm-global/bin/claude /usr/local/bin/claude

# python3, for the same reason and one more caller: the remote frontend's relay
# installer runs node-gyp in a non-login non-interactive shell -- what
# `ssh host '<command>'` gets -- to build a native `node-pty`, because Linux has
# no prebuild. With no python3 reachable there the frontend connects with a
# working file tree and dead terminals, which is a confusing way to lose the one
# feature that matters.
#
# What this symlink actually decides is *which* python3 wins, and that is worth
# being exact about. Debian's python3-minimal already puts 3.11 at
# /usr/bin/python3, which is on sshd's built-in PATH, so the build was never
# failing for want of an interpreter -- measured off a relay's own
# build/config.gypi, which recorded `"python": "/usr/bin/python3"`.
# /usr/local/bin precedes /usr/bin in every default PATH there is, so this hands
# node-gyp the python feature's interpreter instead, which is the one the
# interactive environment already uses. Both were checked against this package:
# `node-gyp configure` succeeds under either, and a full rebuild under the
# feature's 3.14 produces the same pty.node as 3.11 did.
#
# Created before the feature that provides the target: features are layered on
# top of this image, so this is a dangling symlink until the build finishes and
# a working one afterwards. A dangling link is not executable, so a PATH search
# skips it and falls through to /usr/bin/python3 -- which is why getting this
# wrong degrades rather than breaks. orca-target.sh checks over ssh regardless.
RUN ln -sf /usr/local/python/current/bin/python3 /usr/local/bin/python3
USER node


# Copy and set up firewall script
COPY init-firewall.sh /usr/local/bin/
USER root
RUN chmod +x /usr/local/bin/init-firewall.sh && \
  echo "node ALL=(root) NOPASSWD: /usr/local/bin/init-firewall.sh" > /etc/sudoers.d/node-firewall && \
  chmod 0440 /etc/sudoers.d/node-firewall

# The sshd the remote frontend attaches to. Copied in rather than read off
# /workspace: the bind mount is writable by whatever is running in the container,
# and a config that decides who may log in should not be.
COPY sshd-remote.conf /etc/ssh/sshd_remote_config
COPY start-sshd.sh /usr/local/bin/
RUN chmod 0644 /etc/ssh/sshd_remote_config && \
  chmod +x /usr/local/bin/start-sshd.sh && \
  echo "node ALL=(root) NOPASSWD: /usr/local/bin/start-sshd.sh" > /etc/sudoers.d/node-sshd && \
  chmod 0440 /etc/sudoers.d/node-sshd && \
  rm -f /etc/ssh/ssh_host_*

# The rm above is the point of that line, not tidying. The openssh-server
# postinst runs ssh-keygen -A, so a host key would otherwise be baked into the
# image and shared by every container built from it. start-sshd.sh regenerates
# what is missing at first start, which makes it per container instead.

# The volume-ownership repair, on the same sudo arrangement, for a runtime that
# does not seed a named volume from the image. See fix-volume-perms.sh.
COPY fix-volume-perms.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/fix-volume-perms.sh && \
  echo "node ALL=(root) NOPASSWD: /usr/local/bin/fix-volume-perms.sh" > /etc/sudoers.d/node-volperms && \
  chmod 0440 /etc/sudoers.d/node-volperms

# Two things about the node account, and the second is not cosmetic.
#
# sshd reads the login shell out of /etc/passwd and ignores ENV SHELL, so without
# --shell an ssh session lands in bash -- no ~/.zshrc, no gh-auth.sh, no gh.
#
# `usermod -p '*'` unlocks the account. node ships shadow-locked, and OpenSSH's
# allowed_user() refuses a locked account for *public key* auth when UsePAM is
# off, which sshd-remote.conf sets. Every visible thing is correct --
# authorized_keys byte-exact, .ssh 700, the file 600 -- and the client still sees
# a bare `Permission denied (publickey)`. A test fixture that logs in as root
# never meets this, because root is never locked; it appears the moment a normal
# user is the one logging in. `*` means "no password accepted", not "no password
# required": password auth stays off, and PasswordAuthentication no in
# sshd_remote_config says so a second time.
RUN usermod --shell /bin/zsh node && \
  usermod -p '*' node

# /workspace is a virtiofs bind mount and the owner it reports for the mount root
# flaps between node and root from one syscall to the next. git's ownership check
# believes it, so `fatal: detected dubious ownership in repository at
# '/workspace'` arrives intermittently and on some commands and not others --
# measured at 48 of 50 `git status` runs in one phase and 0 of 50 in another.
# Nothing is wrong with the checkout; writes succeed throughout.
#
# --system rather than a session's --global: /home/node is not one of the mounted
# volumes, so a global config is lost on the next rebuild.
#
# `*` and not `/workspace`, which is what this said first and which would have
# left the case that matters broken. git reports and matches the *worktree's* own
# path, so a linked worktree under .claude/worktrees/ is a separate entry --
# measured with `GIT_TEST_ASSUME_DIFFERENT_OWNER=1`, which forces the check:
# `/workspace` alone still refuses the worktree, the exact worktree path is
# accepted, and neither `/workspace/*` nor `/workspace/.claude/worktrees/*`
# matches anything, because git 2.39 interpolates these paths but does not glob
# them -- `*` on its own is the only wildcard the documentation defines. Worktree
# names are made per session, so there is no fixed list to enumerate, and the
# agent workflow runs inside those worktrees rather than in /workspace.
#
# The check exists to stop another user's repository running its hooks as you.
# This container has one user, and the only things mounted are that user's own
# workspace and their own state volume, so there is no second owner to be
# confused with -- which is the trade being made, not an oversight.
RUN git config --system --add safe.directory '*'
USER node
