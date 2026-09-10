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
RUN SNIPPET="export PROMPT_COMMAND='history -a' && export HISTFILE=/commandhistory/.bash_history" \
  && mkdir /commandhistory \
  && touch /commandhistory/.bash_history \
  && chown -R $USERNAME /commandhistory

# Set `DEVCONTAINER` environment variable to help with orientation
ENV DEVCONTAINER=true

# Create workspace and config directories and set permissions.
#
# The memory path is here for a mount this repository does not have yet: a named
# volume is initialised from the image path it covers, ownership included, so a
# path the image does not have comes up root-owned and unwritable by node — with
# no warning, and sudo here is restricted to init-firewall.sh. decoder hit exactly
# that on a shared memory volume: five files failed to restore and the script died
# on a bare Permission denied. Creating it node-owned here means a volume mounted
# over it inherits that instead. /home/node/.config/gh is the same case and does
# have a mount: `gh auth login` writes hosts.yml there, and a root-owned volume
# would fail that write with a bare Permission denied.
RUN mkdir -p /workspace /home/node/.claude/projects/-workspace/memory \
  /home/node/.config/gh && \
  chown -R node:node /workspace /home/node/.claude /home/node/.config

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
  -a "export PROMPT_COMMAND='history -a' && export HISTFILE=/commandhistory/.bash_history" \
  -x

# Install Claude
RUN npm install -g @anthropic-ai/claude-code@${CLAUDE_CODE_VERSION}


# Copy and set up firewall script
COPY init-firewall.sh /usr/local/bin/
USER root
RUN chmod +x /usr/local/bin/init-firewall.sh && \
  echo "node ALL=(root) NOPASSWD: /usr/local/bin/init-firewall.sh" > /etc/sudoers.d/node-firewall && \
  chmod 0440 /etc/sudoers.d/node-firewall

# The sshd `cmux ssh` attaches to. Copied in rather than read off /workspace: the
# bind mount is writable by whatever is running in the container, and a config
# that decides who may log in should not be.
COPY sshd-cmux.conf /etc/ssh/sshd_cmux_config
COPY start-sshd.sh /usr/local/bin/
RUN chmod 0644 /etc/ssh/sshd_cmux_config && \
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
# --shell a cmux ssh session lands in bash -- no ~/.zshrc, no gh-auth.sh, no gh.
#
# `usermod -p '*'` unlocks the account. node ships shadow-locked, and OpenSSH's
# allowed_user() refuses a locked account for *public key* auth when UsePAM is
# off, which this config sets. Every visible thing is correct -- authorized_keys
# byte-exact, .ssh 700, the file 600 -- and the client still sees a bare
# `Permission denied (publickey)`. cmux's own fixture logs in as root, which is
# never locked, so this only appears once you use a normal user. `*` means "no
# password accepted", not "no password required": password auth stays off, and
# PasswordAuthentication no in sshd_cmux_config says so a second time.
RUN usermod --shell /bin/zsh node && \
  usermod -p '*' node
USER node
