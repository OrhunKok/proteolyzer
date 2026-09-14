# The sandbox, without VS Code

Two reasons, and the order matters when weighing a replacement. The first is
that the environment is a file in git, so another machine is a clone rather than
an afternoon. The second is that a firewall stands between Claude Code and
everything that is not GitHub, PyPI, npm or the Anthropic API — see
`init-firewall.sh`, which is why the config asks for `NET_ADMIN`.

Most things that claim to replace this replace the second reason only. See
[Moving to another machine](#moving-to-another-machine).

## Setup

The runtime is Apple `container`; Docker still works and is the fallback path.

```bash
brew install wcgomes/tap/adevcontainer        # once — needs Apple container
adevcontainer doctor                          # checks the runtime is usable

./.devcontainer/build.sh                      # build the image
./.devcontainer/orca-target.sh                # prepare it, print the ssh target
```

That last one prints the four fields to type into the frontend's Add Target form,
once. Day to day:

```bash
./.devcontainer/up.sh                # a zsh in the container
./.devcontainer/up.sh claude         # straight into Claude Code
./.devcontainer/orca-target.sh       # re-check the target after a rebuild
REBUILD=1 ./.devcontainer/up.sh      # replace the container
./.devcontainer/build.sh             # after any Dockerfile change
```

`build.sh` is separate because `devcontainer.json` says `image:` rather than
`build:` — the one change that makes this run on Apple `container` at all, since
`adevcontainer` rejects Dockerfile builds. The cost is that a Dockerfile edit no
longer takes effect on its own; build first. The gain is that the same config
serves both runtimes, because the devcontainer CLI reads `image:` too.

`up.sh` picks `adevcontainer` when it is installed and the devcontainer CLI
otherwise. The ssh path is Apple-only, because that is where the container has
its own address and the whole published-port dance disappears.

## Names, and giving the container one

`name` in `devcontainer.json` is the project — here `proteolyzer` — not a
description of what the container is for. adevcontainer turns it into both the
create name and the DNS hostname, so `"Claude Code Sandbox"` gave every project a
container called `claude-code-sandbox` — indistinguishable from the next one in
Orchard or `container list`, and useless as a hostname.

It is `${localWorkspaceFolderBasename}` rather than a literal, so the file is
identical in every project and the folder name is the only thing that decides.
See [Migrating another project onto this](#migrating-another-project-onto-this).

A `buildkit` container appearing beside it is Apple's own builder for
`container build`, not one of these. It comes and goes.

**A stable name matters** beyond looks: every rebuild gets a fresh address, so a
saved ssh target pinned to an IP goes stale the next time you rebuild. This project is `proteolyzer.adevcontainers.local`,
and it resolves — but it did not for a while, and the reason is worth keeping.

### Why the record was missing, and what a working one costs

It was missing, established on 2026-09-10 by control:

```bash
container run -d --rm --name dnstest docker.io/library/alpine sleep 300
dig @127.0.0.1 -p 2053 +short dnstest.adevcontainers.local   # 192.168.64.3
dig @127.0.0.1 -p 2053 +short proteolyzer.adevcontainers.local   # nothing
```

Everything else was ruled out first: the domain in `container system dns list`
and in `container system property list`, a correct
`/etc/resolver/containerization.adevcontainers.local` (the filename prefix is
cosmetic — macOS reads the `domain` directive inside), the service restarted, and
the container recreated afterwards and running as `proteolyzer`. `dig` came back
empty rather than refused, so the service was reachable and simply had no record.

That was first read as *adevcontainer's containers do not get DNS records, while
`container run --name` ones do*. *Apple's source does not support that reading*,
and it is worth not acting on: it makes the ssh alias below a workaround for
someone else's bug, which is the kind of thing that never gets revisited.
`container create` and `container run` build their configuration through one
shared helper, and the name a container registers under is decided there, once:

```swift
// apple/container — Sources/Services/ContainerAPIService/Client/Utility.swift
config.networks = try getAttachmentConfigurations(   // :208
    containerId: config.id, …, dnsDomain: containerSystemConfig.dns.domain)
// and inside it: fqdn = "\(containerId).\(dnsDomain)." — or nil, if that domain is nil
```

With that domain nil the container attaches as the bare `proteolyzer`, and a bare
hostname is what the resolver then has to answer with. **The value is read once,
at create time, and baked into the stored config** — so setting the domain
afterwards changes nothing until the container is *recreated*, not merely
restarted. Create-versus-run is not the variable; created-before-versus-after the
domain was set is.

`container system dns create <domain>` does not set it. That writes
`/etc/resolver/…` and the pf rule, which points macOS *at* the container DNS
service and stops there. What containers register under is the `dns.domain`
system property:

```toml
# ~/.config/container/config.toml
[dns]
domain = "adevcontainers.local"
```

**That file is not the one anything reads.** It is a source, copied to
`~/Library/Application Support/com.apple.container/config/config.toml`, and
`container create` reads only the copy (`Application.swift:149-158`). The copy is
made by exactly one caller — `container system start`
(`SystemStart.swift:76` → `ConfigurationLoader.copyConfigurationToReadOnly`). So
editing the user config does nothing at all until the service is restarted, and
then does nothing to containers that already exist. Two steps, both easy to miss,
and neither reports anything.

`container system property list` reads the copy, so it answers "is this live?"
rather than "did I write it?" — which is why it is the right thing to check, and
why checking it on its own is not enough to say when it *became* live.

**One line reads the answer off a running container**, which is what the
investigation above was missing:

```bash
grep domain /etc/resolv.conf   # nothing == no domain when this container was created
```

That line comes from the same value (`RuntimeService.swift:1198`, via
`DNS.resolvConf` in containerization, which emits `domain …` whenever the domain
is non-nil). `/etc/hosts` cannot answer it — the guest hostname is deliberately
truncated to its first label (`:1190`), so it reads `proteolyzer` either way.

On this container, on 2026-09-10, there is no `domain` line, while
`container system property list` on the Mac shows `domain = "adevcontainers.local"`
the same day. Both are true and they are not in conflict: **the property is live
now and was not when this container was created.** Nothing rewrites the file
afterwards — no script here touches it, it is a plain file rather than a mount,
and its mtime predates the boot — so it still reads as it was written at create.

The reading is not ambiguous either. `config.dns` was present with empty
nameservers, which is the branch that fills them from the attachment gateway
(`RuntimeService.swift:229-238`) — hence `nameserver 192.168.64.1`. Had the
domain been set, the same struct would have carried it. `--no-dns` would have
produced neither line.

**So the fix was a rebuild, and nothing else** — confirmed on 2026-09-11. The
container was recreated with the property already live, and this time it came up
with the domain and a record to match:

```bash
# inside the new container
cat /etc/resolv.conf                                       # nameserver 192.168.64.1
                                                           # domain adevcontainers.local
dig @192.168.64.1 +short proteolyzer.adevcontainers.local   # 192.168.64.6
```

Nothing was changed to achieve that. No flag, no `runArgs`, no config edit — the
same image and the same `devcontainer.json`, created once more after the property
had taken effect. Which is the whole claim: the name was never adevcontainer's to
give or withhold.

So `container system dns create` *and* `dns.domain` *and* a rebuild, in that
order, is what a working name costs — and the middle one needs
`container system start` before it counts.

**Nothing is written to `~/.ssh/config`.** The attach script used to add a `Host`
block there per project, and that was wrong twice over. A file every project
appends to is shared state, and one container has no business knowing which
others exist — the same objection that collapsed three volumes into one.

It also did not buy what it was sold on. The claim was that a launch command
pinned to a pane would go stale when a rebuild changed the address. It would not:
that command runs the script, which recomputes the address every time, so nothing
pinned ever held an IP. And the block's `HostName` was static until the script
rewrote it, so anything that reconnected on its own was exactly as stale as a
bare address would have been. It bought a prettier string in a log.

What runs now: the name when this Mac can resolve it, the address when it
cannot, passed on the command line either way.

```bash
dscacheutil -q host -a name proteolyzer.adevcontainers.local
```

`dscacheutil` and not `dig @127.0.0.1 -p 2053` — the second proves the record
exists, the first proves *this Mac* will find it, and only that decides whether
an attach works. `DEVCONTAINER_SSH_HOST` names the target outright,
`CONTAINER_DNS_DOMAIN` just its suffix.

A `ProxyCommand` resolving the address at connect time was tried before either
and is tidier in principle. A frontend that bootstraps a daemon on the remote
could not get through it — `failed to query remote platform: Connection closed by
UNKNOWN port 65535`, where `UNKNOWN` is ssh not knowing its peer because a proxy
is in the way. That bootstrap does much more than open a session: platform probe,
binary upload, forward.

### What is idiomatic here and what is a workaround

Worth separating, because the workarounds are the parts to delete the day
upstream fixes them.

**As Apple intends it.** The image is ordinary OCI from `container build`.
`capAdd: ["NET_ADMIN"]` is the documented spelling — Apple's
`runtime-configuration.md` gives `container run --cap-add NET_ADMIN` verbatim.
State is a named volume through `container volume`. And the container is reached
on **its own address with nothing published**, which is Apple's networking model
rather than a Docker habit carried over; the published-port version of this file
was the less idiomatic one.

**Workarounds, and what each is for.** `fix-volume-perms.sh` exists because Apple
does not seed a named volume from the image path it covers
([#729](https://github.com/apple/container/issues/729)). `build.sh` deletes
derived images because adevcontainer's derived tag hashes the config and not the
base, so a rebuilt base is silently ignored. Those two are gaps in a specific
tool, not disagreements with Apple's design, and each should be removed rather
than maintained once it closes.

The ssh alias is gone: it was neither a gap in anything nor worth the
shared file it lived in. See above.

**One thing left that is more workaround than it needs to be.**
`sshd-remote.conf` sets `UsePAM no`, which is why the `node` account has to be
unlocked with `usermod -p '*'` — OpenSSH refuses a locked account for public key
auth when PAM is off. `UsePAM yes` is Debian's own default and would need
neither. It is left alone because the current arrangement is proven and swapping
an authentication path that works for one that is merely more standard is a poor
trade; but that is the reason, not an argument that this is better.

The VS Code extension still works — `customizations.vscode` is read when it
attaches — but nothing depends on it.

## More than one project at once

`adevcontainer exec` opens an interactive **"Select a container:"** picker when
more than one managed container is running. A script that execs four times
therefore stops four times, and acts on whichever row happens to be highlighted
— not on the project it was started from.

That is not theoretical. With `pinpoint` and `proteolyzer` both up, an attach
from `proteolyzer` wrote the ssh key into one container, read sshd's pid from
the other, and reported the wrong address as the project's, so the alias for
`proteolyzer` pointed at `pinpoint`'s IP and Claude opened with the wrong
history. The hangs were the pickers waiting.

So every `exec` in `up.sh` and `ssh-target.sh` passes `--name`, taken from the
`containerId:` line of the tool's own output rather than guessed from the
directory. The host that gets reported follows the same value, so it names the
container it actually reaches even if `name` in `devcontainer.json` and the
folder ever disagree.

Worth knowing when you run these by hand too: `adevcontainer stop`,
`adevcontainer exec` and friends will all prompt. Pass `--name <project>`.

## Skipping permission prompts

The palette entries run `claude --dangerously-skip-permissions`, which is what
the firewall is for: a sandbox that cannot reach anything but GitHub, PyPI, npm
and the Anthropic API is the environment where bypassing prompts is a reasonable
trade rather than a reckless one. Anthropic's own reference devcontainer exists
for the same purpose.

The flag rather than the setting, deliberately. `permissions.defaultMode` set to
`bypassPermissions` does the same thing, but since v2.1.257 it is **only honoured
from user or managed settings** — a `.claude/settings.json` committed to a
repository is read and then ignored for this key, silently. If you want typed
`claude` sessions to bypass as well, it has to go in the container's *user*
settings, which is inside the volume and therefore survives rebuilds and travels
with `state.sh`:

```bash
# inside the container
jq '.permissions.defaultMode = "bypassPermissions"' \
   "$CLAUDE_CONFIG_DIR/settings.json" > /tmp/s && mv /tmp/s "$CLAUDE_CONFIG_DIR/settings.json"
```

Two things this does not do. It is not protection against prompt injection —
nothing here is. And the firewall bounds *where* the agent can reach, not what it
can do with what is already inside: the `gh` credential in the volume can push to
your repositories. Deny rules still apply in every mode, `bypassPermissions`
included, so they are the place to put anything that should stay impossible.

## What VS Code was also doing, silently

Two things came from the editor rather than from `devcontainer.json`, and both
have to be replaced explicitly when the container is started from a terminal.

**The GitHub credential.** VS Code forwards a credential from the host by
installing a git credential helper inside the container that calls back to the
editor over an RPC pipe. That is why `git push` used to work with no setup, and
`gh-auth.sh` existed only to lend the same credential to `gh`, which does not
read git helpers. `devcontainer up` has no editor and forwards nothing: the pipe
is absent, `git credential fill` returns empty, and *both* tools are dead — `gh`
loudly, `git push` with whatever error the stale helper produces.

So the container holds its own credential now. Once per volume:

```bash
gh auth login          # inside the container; device flow, github.com is allowed
```

That login lands in `~/.state/gh`, inside the one named volume keyed on the
directory name, so it survives a rebuild the way the history and the Claude
config do.
`gh-auth.sh` notices it on the next shell and runs `gh auth setup-git`, which
points git at the same credential — one login, both tools, no host involved. The
old forwarding path stays as the fallback for a VS Code session.

One edge: `gh auth setup-git` is skipped when a credential helper is already
configured, because two helpers is a coin toss over which answers. If you attach
VS Code to a container and then go back to the terminal, its helper is left
behind pointing at a pipe that is gone. `git config --global --unset-all
credential.helper` clears it, or rebuild.

**Formatting on save.** `editor.formatOnSave` and the eslint fixer in
`customizations.vscode` do nothing without the editor. `make lint` and
`.pre-commit-config.yaml` are what enforce formatting here, and they run in the
container either way.

**Where you edit is not one of them**, though it looks like it should be. The
Dev Containers extension puts the editor *inside* the container, so losing it
reads as losing somewhere to write code. But `workspaceMount` is a bind: the
files are on the Mac the whole time,
and `/workspace` is a view of them. Any native editor opens the checkout
directly, with no container in the path. What belongs in here is running things
— the agent, the tests, anything that should be behind the firewall — not
typing. That is a better arrangement than the one being replaced, not a
casualty of it.

Three entries in `init-firewall.sh` — `marketplace.visualstudio.com`,
`vscode.blob.core.windows.net`, `update.code.visualstudio.com` — exist so the VS
Code server and its extensions can install themselves inside the container. They
are dead weight once nothing attaches. They are left in because they cost
nothing and are the difference between VS Code working and VS Code hanging on
the day you want it back.

## The remote frontend

The container runs sshd, and an app on the Mac treats it as a remote host. That
is the whole arrangement: **a frontend needs a host, a port and an identity file,
and `ssh-target.sh` produces those for whoever asks.** `orca-target.sh` is a thin
thing on top of it, which is why trying a different app costs an afternoon rather
than a migration.

That leaves two ways in, and they are not equivalent:

| | `up.sh` (`devcontainer exec`) | `orca-target.sh` (ssh) |
|---|---|---|
| extra surface in the image | none | `openssh-server` on port 2222 |
| runtime | Apple `container` or Docker | Apple `container` only |
| file tree, editor, diff view | no | yes, on the container's filesystem |
| terminal survives the app quitting | no | yes, leased by a relay on this side |
| sftp file drop | no | yes |
| `TERM` in the container | pinned to `xterm-256color` | the host's own, `COLORTERM` included |

Use `up.sh` for a shell. Use the ssh path for the workspace you actually work in.

### How the ssh path works

What this repository adds for it:

- `openssh-server` in the image, and `sshd-remote.conf` — port 2222, no root
  login. `AllowTcpForwarding yes` is load-bearing, not boilerplate: a frontend
  that leases terminals through a relay on this side reaches that relay over the
  ssh connection, and asks with `ExitOnForwardFailure=yes`, so refusing it fails
  the connection rather than degrading it.
- `start-sshd.sh`, reachable through `sudo` for the same reason
  `init-firewall.sh` is, and run from `postStartCommand` after it — the firewall
  flushes the tables, so anything holding a connection wants to start after it.
  Host keys are generated at first start rather than baked into the image, so two
  containers from one image do not share one.
- `usermod --shell /bin/zsh node`. sshd reads the login shell out of
  `/etc/passwd` and ignores `ENV SHELL`, so without this an ssh session lands in
  bash — no `~/.zshrc`, so no `gh-auth.sh`, so no `gh`. Everything else names
  zsh explicitly and is unaffected.
- `usermod -p '*' node`, because `UsePAM no` makes OpenSSH refuse a locked
  account for public key auth. See the workaround note above.

**No port is published.** Under Apple `container` every container has its own
address, reachable from the host, so there is nothing to discover and nothing
bound to loopback — which is also why this path is Apple-only. `capAdd` replaced
`runArgs` once adevcontainer turned out to reject `-p` entries outright, and the
ssh path stopped needing a published port in the same move.

Starting sshd unconditionally is safe: nothing is published, password auth is
off, and no `authorized_keys` exists until `ssh-target.sh` writes one. Until you
ask for it, it listens and refuses everything.

The key is a dedicated one, `~/.ssh/devcontainer`, not the one that talks to
GitHub — giving a hop into a sandbox a key with any other reach is how a sandbox
stops being one. Agent forwarding is explicitly off for the same reason. Host key
checking is off because a host key is generated per container and the address
changes on every rebuild, so `known_hosts` could only ever reject a container it
had seen before at an address something else now answers on. What bounds this is
the address being on Apple's container network and the key being one the script
made.

### A second app on the same sshd

The Claude Code desktop app has an SSH environment — its documentation names dev
containers as a target — and it asks for exactly what `ssh-target.sh` already
produces: a host, a port, and an identity file.
`node@proteolyzer.adevcontainers.local`, `2222`, and `~/.ssh/devcontainer`. It
installs Claude Code on the remote itself and uses the remote
`/home/node/.claude`, which is the volume.

Two things to expect if you try it. There is no terminal panel in a remote
session, so it wants a terminal beside it rather than replacing one. And SSH
mode sets `CLAUDE_CODE_PROVIDER_MANAGED_BY_HOST=1` with its own
`ANTHROPIC_BASE_URL` and token, overriding provider configuration on the remote
— which costs nothing here, and would matter to a container pointed at Bedrock
or Vertex.

### Orca

[Orca](https://github.com/stablyai/orca) is an MIT-licensed Electron app from
Stably AI — parallel agents, each in its own worktree — with an SSH mode carrying
a file tree, a fuzzy finder, an editor and a diff view that all operate on the
remote filesystem. `orca-target.sh` prepares the container and prints the four
fields its Add Target form wants.

The reason to care is not the editor. It is this, from its `docs/site/content/docs/ssh.mdx`:

> Closing the desktop app no longer kills your remote PTY sessions. Remote
> terminal sessions are leased through the relay running on the remote host, so
> they survive Orca closing on your laptop. When you reopen the app and reconnect
> to the target, leased PTYs are restored to their tabs in the **attached**
> state, with their scrollback intact.

The session keeps running on the container instead of being relaunched from a
command the host remembered, which is the right end to solve it from. "Keep
terminals alive until reset" is on by default for every target.

**What it needs from the image**, and it is a hard requirement: the relay builds
a native `node-pty` on the remote and Linux has no prebuild — the package ships
`prebuilds/` for darwin and win32 only — so the container needs `make`, a C++
compiler and `python3`. Without them Orca still connects — files, git and the
editor all work — and remote terminals simply do not, which is a confusing way to
lose the only feature that matters. `node:20` brings the first two from
`buildpack-deps`. `python3` is present twice over: Debian's `python3-minimal` at
`/usr/bin/python3`, and the python feature's, which installs outside
`/usr/local/bin` and is therefore invisible to the non-login shell Orca's
installer runs in — hence the symlink in the Dockerfile, next to the one `claude`
needs for the same reason. `orca-target.sh` checks all three over ssh rather than
trusting any of this.

**Two things to set.** Point the target at the DNS name, not the address: an Orca
target is a saved host entry, and the address changes on every rebuild. And set
the repo's worktree base path to something **relative** — `.claude/worktrees` is
what this repository already uses. Orca ignores an absolute base path for an ssh
repo and falls back to a sibling of the repo (`join(repoPath, '..', …)` in
`src/main/ipc/worktree-logic.ts`); the repo here is `/workspace`, so a sibling is
the container's own root filesystem, which is invisible from the Mac, outside
the state volume, and destroyed by the next rebuild.

**It rewrites `~/.claude/settings.json`.** A first connect installs Orca's own
hook block there and keeps nothing else — `model`, and any hook already in the
file, are gone afterwards, silently. Copy the file aside first. An older copy is
recoverable from a state tarball:

```bash
tar -xzOf devcontainer-state.tar.gz --wildcards '*/.claude/settings.json'
```

One more that has not been chased: the relay keeps terminal history under
`~/.orca-remote/` in the container's home, which is *not* the state volume, so
scrollback does not survive a rebuild the way the Claude and `gh` state does.
The fix is the same shape as the rest — a path under `/home/node/.state` — if it
turns out to matter.

### Notifying from inside the container

An OSC escape sequence needs nothing installed:

```json
{
  "hooks": {
    "Notification": [
      {
        "hooks": [
          {
            "type": "command",
            "command": "printf '\\033]9;Claude needs input\\007' > /dev/tty"
          }
        ]
      }
    ]
  }
}
```

`> /dev/tty` is the part that matters — hook stdout is captured by Claude Code,
so a sequence merely printed never reaches the terminal. That settings.json
belongs in `/home/node/.claude`, which is a persisted volume — and which the note
above about Orca rewriting it applies to.

### Other answers to the same problem

[ccmanager](https://github.com/kbwo/ccmanager) runs the agent session inside the
devcontainer as a first-class feature, with the manager on the host. The trade is
a session manager instead of a terminal — no browser pane, no splits, no socket
API.

### A different runtime under all of this

Docker is not the only way to get a Linux container on a Mac, and the one that
suits this best is Apple's own — every container gets an address reachable from
the host, which deletes the published port and the port-discovery step from the
ssh path rather than adding to them. `NET_ADMIN` is supported, so the firewall
survives the move. What does not work yet is `build.dockerfile`, which is what
this repository uses. [APPLE.md](./APPLE.md) has the evidence, the config to swap
in, and the one upstream issue to watch. OrbStack is the meanwhile option: it
makes Docker faster without making it different, so nothing here changes.

### The workspace mount lies about who owns it

`/workspace` is a virtiofs bind mount, and the owner it reports for the mount
root flaps between `node` and `root` from one syscall to the next. git's
ownership check believes it:

```
fatal: detected dubious ownership in repository at '/workspace'
```

Intermittently, and on some commands and not others, which is what makes it read
like anything but a mount problem. `git log` succeeds while `git fetch` in the
same second fails — a fetch runs `git rev-list` and `git maintenance` as
subprocesses and each re-runs the check, so it has more chances to land on a bad
sample. Measured at 48 of 50 `git status` runs failing in one phase and 0 of 50
in another. Nothing is wrong with the checkout: writes succeed throughout,
including while it reads `root:root 700`.

The Dockerfile answers it with `git config --system --add safe.directory '*'`.
`--system` rather than a session's `--global`, because `/home/node` is not one of
the mounted volumes and a global config is lost on the next rebuild.

**`*` rather than `/workspace`**, which is what this said first and which would
have left the case that matters broken. git matches — and reports — the
*worktree's* own path, so each linked worktree under `.claude/worktrees/` is a
separate entry. Forcing the check with `GIT_TEST_ASSUME_DIFFERENT_OWNER=1` shows
it plainly:

```
$ GIT_TEST_ASSUME_DIFFERENT_OWNER=1 git -c safe.directory=/workspace status
fatal: detected dubious ownership in repository at
'/workspace/.claude/worktrees/salvage-111'
```

The exact worktree path is accepted; `/workspace/*` and
`/workspace/.claude/worktrees/*` match nothing, because git 2.39 interpolates
these paths but does not glob them — `*` alone is the only wildcard the
documentation defines. Worktree names are made per session, so there is no fixed
list to enumerate, and the agent workflow runs *inside* those worktrees rather
than in `/workspace`: an entry for `/workspace` only would have fixed the case
nobody was hitting.

The check exists to stop another user's repository running its hooks as you. This
container has one user, and the only things mounted are that user's own workspace
and their own state volume, so there is no second owner to be confused with. That
is the trade, and it is worth stating rather than implying.

## Migrating another project onto this

The procedure, in order, for a project still on Docker and the VS Code extension.
It has been run three times; every step below exists because one of them went
wrong.

**Before anything, look at what you are about to delete.**

```bash
cd ../otherproject
git status                    # commit or stash first -- step 2 is destructive
ls .devcontainer              # anything here that is not the stock template?
```

The stock Claude template is replaced wholesale by this. A project that added
packages to its own `Dockerfile` wants those carried over by hand, and this is
the only moment you will remember to look.

**1. Start Docker Desktop.** The project's history lives in its old Docker
volumes, and step 4 reads them. Skip this and the container comes up with no
sessions — which is recoverable, but only if you notice.

**2. Replace the config.**

```bash
rm -rf .devcontainer
cp -R ../thisproject/.devcontainer .
```

The `rm -rf` is not tidiness, it is the whole trap. `cp -R src dest` copies
*into* `dest` when `dest` exists, so a project that already has a
`.devcontainer/` gets `.devcontainer/.devcontainer/` and keeps its old
`devcontainer.json`. Everything then proceeds as though the copy worked, until
adevcontainer reads the old file and says `Dockerfile build is not supported` —
which reads as a fault in the new setup and is the old config still sitting
there. `cp -R ../thisproject/.devcontainer/. .devcontainer/` is the alternative,
and for a project under git the delete is safe anyway: `git checkout
.devcontainer` brings the originals back.

Nothing in that directory names a project, so there is nothing to edit. Check
rather than trust:

```bash
grep -rn thisproject .devcontainer          # expect nothing
```

**3. Build, once per machine rather than once per project.**

```bash
./.devcontainer/build.sh
```

Every project runs `claude-devcontainer:local`, so after the first this is a
no-op and takes seconds.

**4. Start it, and watch what it says.**

```bash
./.devcontainer/up.sh
```

Plain `up.sh`, not `REBUILD=1`: there is no container for this folder yet, and
nothing to rebuild. You are looking for

```
state.sh: adopting otherproject's earlier state from docker
state.sh: adopted. Anything already in the volume was left alone.
```

`nothing to adopt` means Docker was reachable and the old volumes are not there
— check `docker volume ls | grep claude-code`. A complaint that Docker is not
running means step 1 was skipped; start it and run `up.sh` again, which is safe
because the volume is only marked done when the answer was knowable.

**5. Log in, once per project**, because each has its own volume:

```bash
gh auth login
claude                       # paste the URL with ⌘V; selecting it truncates it
```

**6. Confirm the history actually arrived**, rather than assuming:

```bash
# on the Mac
docker run --rm -v claude-code-config-otherproject:/v alpine \
  sh -c 'ls /v/projects/-workspace/*.jsonl 2>/dev/null | wc -l'
# inside the container
ls ~/.state/claude/projects/-workspace/*.jsonl | wc -l
```

Those are the conversation transcripts and the numbers should match. `-workspace`
is the right key for every project: Claude records a project by the path it ran
at, and that is `/workspace` in every container.

**7. Commit it**, or a stray `git checkout` takes the whole setup with it.

Afterwards, run `./.devcontainer/orca-target.sh` and add the target it prints.
The old Docker volumes can be deleted once you trust the copy — not before, since
until adopt has run they are the only one.

## Bringing a project's earlier state with it

A project that has been used with Claude before copying this in has history
worth keeping, and it is not in the new volume — it is in whatever the project
used previously, usually the three `claude-code-*` Docker volumes from the stock
template. Nothing looked for it, so such a project opened with no sessions and
the only remedy was a hand-run sequence. That was a gap in this setup rather
than a thing to do per project.

`up.sh` now runs `state.sh adopt` once per volume. It finds the legacy volumes in
whichever runtime has them — the same one, or Docker when the live container is
Apple's — and copies in **only what the volume does not already have**, so a
credential from a login you just completed is never overwritten. `cp -n` and
`tar --skip-old-files` are both missing from busybox, so that no-clobber is
written out explicitly.

It is idempotent by marker: `/.adopted` in the volume, written whether or not
anything was found, so a project with nothing to adopt does not pay for the
search on every start. Running `./.devcontainer/state.sh adopt` by hand does the
same thing, and the cross-runtime path needs Docker running to read the old
volumes.

## Moving to another machine

Most of it already travels, and it is worth being exact about which part does
not:

| | travels | how |
|---|---|---|
| the environment itself | yes | `Dockerfile` and `devcontainer.json`, in git |
| the workspace | yes | a bind mount, so wherever you cloned |
| the ssh key and host key | yes | `ssh-target.sh` and `start-sshd.sh` make them if absent |
| **the three named volumes** | **no** | `state.sh` |

The volumes are the gap, and one of them matters: `/home/node/.claude` holds
settings, project state and the agent's memory directory. `state.sh export`
writes them to a tarball you copy across and `state.sh import` puts them back.
The `gh` volume is excluded unless you ask for it, because including it writes a
GitHub token into that tarball in plaintext and `gh auth login` is one command
on the far side. Volume names are keyed on the directory basename, so the
checkout has to be named the same over there — the same constraint
`devcontainer.json` already documents.

**The macOS-only parts are the frontend, not the environment.** Neither the
desktop app nor Apple `container` runs on Linux — but neither of them is the
environment. The image is OCI, so Docker, podman, containerd and Apple
`container` all take it and none of them cares which one built it; the runtime
underneath is interchangeable and picking a Mac-native one costs nothing here.
What is runtime-specific is the two ways in, and that is why there are two:
`orca-target.sh` is the nice thing on a Mac, `up.sh` is the one that still works
on a Linux box, and neither is the source of truth.

**`publish.sh` is the step that makes migration a pull instead of a build.** It
cross-builds `linux/amd64` and `linux/arm64` with buildx and pushes both to
GHCR under a tag derived from the contents of this directory, so the tag only
moves when the thing that builds the image does. A machine that has the tag
needs neither the Dockerfile nor a build.

That also makes the Apple `container` route cheaper than [APPLE.md](./APPLE.md)
implies, and for a reason worth spelling out: the one change `adevcontainer`
forces — `image:` instead of `build:`, because it hard-rejects Dockerfile builds
— is the same change that makes a pull possible. Those two goals converge rather
than compete. `features` survives the switch untouched; both the devcontainer
CLI and `adevcontainer` derive an image from base plus features at create time.

**`devcontainer.json` stays the source of truth for the same reason.** Every
alternative weighed here trades a standard several tools implement — the CLI,
VS Code, Codespaces, JetBrains, `adevcontainer` — for a bespoke config only its
own tool reads: AgentBox has `agentbox.yaml`, `claude-contained` has a launcher
script, Sculptor has an app. Each is a better product in some direction and all
of them cost the property that made this worth building.

Worth naming because it is a tempting swap: Claude Code's own sandboxing —
macOS Seatbelt plus a domain-allowlisting proxy — is a real answer to the
*second* reason and no answer at all to the first. A Seatbelt policy does not
carry a Python version, a toolchain or an installed package to a new machine.
Good security, not portability.

## Other repositories

`streamlit-DO-MS` and `decoder` carry their own copy of this container, and so
does every other project with a Claude sandbox. Nothing here reaches into them:
a change is applied in each repository by that repository, which is the rule
everywhere else on the account.

Copying is fine for the config and wasteful for the image. **Nothing in this
`Dockerfile` is specific to this repository** — no Python, no project paths, no
`src/`; it is node, the tools, `gh`, Claude Code, the firewall and sshd. The
Python comes from a `features` entry in `devcontainer.json`, which is where the
per-project part belongs. So the same image serves every project, and the honest
description of the status quo is that one generic image is being maintained N
times by hand.

`publish.sh` with a neutral name is the fix, when it is time for one:

```bash
IMAGE=ghcr.io/orhunkok/claude-devcontainer ./.devcontainer/publish.sh
```

Each project then keeps a thin `devcontainer.json` pinning a tag, with its own
`features`, `mounts` and lifecycle. Changing the firewall becomes one edit, one
publish and N pin bumps, instead of N edits that drift.

**Sharing a base does not cost the option to diverge**, which is the usual
objection and the wrong way round. Divergence has three homes, and only the last
needs a Dockerfile: a `features` entry, which is already where this repository's
Python comes from; the `mounts`, `containerEnv` and lifecycle keys; and failing
those, a project Dockerfile that is `FROM` the shared tag plus whatever that
project needs. A project can leave entirely and still inherit the firewall.

**When to do it is a different question from whether**, and the answer is not
yet. Before this branch, this Dockerfile changed twice in three weeks. While it
is moving, N copies are cheaper than a publish step, because an edit is an edit
rather than an edit, a build, a push and N pin bumps.

The signal to switch is making the same edit twice in two repositories, and it
is worth knowing that this has already happened here. `8f5fc5d2`, "Create the
memory path in the image, node-owned", is a fix for a bug **decoder** hit — five
files failing to restore behind a bare `Permission denied` — applied here
preemptively for a mount this repository does not even have. That knowledge
crossed repositories because a person carried it. [CLAUDE.md](../CLAUDE.md)
records the same shape with `CoordinatesMapping`, which existed twice and had
the same imaging-channel bug fixed in both in the same week without either
knowing.

So the changes this file actually attracts are cross-cutting fixes rather than
project-specific features, and that is the profile that eventually wants one
image. A file that changed weekly with per-project tweaks would not.

**What pinning does not answer**, and the reason the bar for switching should
stay high. A pinned tag makes the *bytes* independent — a project publishing a
new tag cannot move a project that did not bump — but the coupling that remains
is in the process:

- Someone has to publish. A Dockerfile change stops being an edit and becomes an
  edit, a build, a push and N pin bumps.
- **It does not reduce the number of pull requests.** A fix all N projects need
  is still N bumps; what is saved is writing the fix N times, not landing it N
  times. So it pays when a fix is hard to write, or when a project would
  otherwise be forgotten entirely — and barely at all when it is a one-line
  change to an obvious place.
- Divergent needs have to be negotiated into one Dockerfile or pushed into a
  `FROM` layer, where today they are simply unrelated.
- A project that goes `image:`-only cannot build from source any more, so the
  registry becomes load-bearing for starting work at all. Keeping the Dockerfile
  as a fallback means maintaining both.

Copies trade a silent-drift risk for none of that. Three projects and a 4 KB
file is a small enough drift risk to keep paying, and the honest summary is that
the two options are closer than the section above makes them sound.

**Pin the content tag, never `latest`**, and [DECISIONS.md](../DECISIONS.md) is
the reason rather than taste. What was torn out on 2026-08-23 was shared
*content that had to be kept current and misled silently when it was not*, and
the test it leaves behind is: *can it be out of date?* A tag like
`:4a0b18e1ff50` cannot — it is immutable, the consumer owns its pin, and a
project on an old one is behind rather than broken. That is the same shape as
two repositories pinning a wheel from here. `:latest` fails the same test for
the same reason: it changes under a project that did not ask it to.
