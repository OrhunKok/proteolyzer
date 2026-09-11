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
brew install --cask cmux                      # once
brew install wcgomes/tap/adevcontainer        # once — needs Apple container
adevcontainer doctor                          # checks the runtime is usable

./.devcontainer/build.sh                      # build the image
./.devcontainer/cmux-attach.sh claude         # a cmux workspace running Claude
```

Day to day:

```bash
./.devcontainer/up.sh                # a zsh in the container
./.devcontainer/up.sh claude         # straight into Claude Code
./.devcontainer/cmux-attach.sh       # the full cmux workspace
REBUILD=1 ./.devcontainer/up.sh      # replace the container
./.devcontainer/build.sh             # after any Dockerfile change
```

`build.sh` is separate because `devcontainer.json` says `image:` rather than
`build:` — the one change that makes this run on Apple `container` at all, since
`adevcontainer` rejects Dockerfile builds. The cost is that a Dockerfile edit no
longer takes effect on its own; build first. The gain is that the same config
serves both runtimes, because the devcontainer CLI reads `image:` too.

`up.sh` picks `adevcontainer` when it is installed and the devcontainer CLI
otherwise. `cmux-attach.sh` is Apple-only, because that is where the container
has its own address and the whole published-port dance disappears.

## Names, and giving the container one

`name` in `devcontainer.json` is the project — here `proteolyzer` — not a
description of what the container is for. adevcontainer turns it into both the
create name and the DNS hostname, so `"Claude Code Sandbox"` gave every project a
container called `claude-code-sandbox` — indistinguishable from the next one in
Orchard or `container list`, and useless as a hostname.

It is `${localWorkspaceFolderBasename}` rather than a literal, so the file is
identical in every project and the folder name is the only thing that decides.
See [Copying this into another project](#copying-this-into-another-project).

A `buildkit` container appearing beside it is Apple's own builder for
`container build`, not one of these. It comes and goes.

**A stable name matters** beyond looks: every rebuild gets a fresh address, so a
cmux workspace or a `cmux surface resume set` command pinned to an IP goes stale
the next time you rebuild. This project is `proteolyzer.adevcontainers.local`,
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

**The ssh alias predates all of that, and still runs.** `cmux-attach.sh` writes a
`Host` block to `~/.ssh/config` with the container's current address, and
rewrites it on every run. It was built when the name did not resolve, and it is
kept because it does not depend on whether the name resolves — the workaround
list below says why that is worth something rather than just redundant.

A `ProxyCommand` resolving the address at connect time was tried first and is
tidier in principle. cmux could not bootstrap its remote daemon through it —
`failed to query remote platform: Connection closed by UNKNOWN port 65535`,
where `UNKNOWN` is ssh not knowing its peer because a proxy is in the way. That
bootstrap does much more than open a session: platform probe, binary upload,
reverse forward. A plain connection to the address was already proven to work, so
the alias stays and the proxy went. The only staleness window is a container
restarted without running this script — and the script is what a
`cmux surface resume set` command runs, so it closes itself.

The block is prepended rather than appended, deliberately: ssh takes the *first*
value it sees for each keyword, so a `Host *` earlier in the file would win on
`IdentityFile` and the right key would never be offered. It is marked with
`# BEGIN cmux-devcontainer …` and replaced in place on each run, so it refreshes
rather than accumulating. Deleting the block opts out.

`CMUX_DEVCONTAINER_HOST` sets the alias, `CONTAINER_DNS_DOMAIN` just its suffix.
Neither needs DNS to be configured at all now.

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

The ssh alias is the odd one out: it is not a gap in anything. A container
registers under whatever `dns.domain` was live when it was created, and for this
one that was nothing. The alias still earns its place — it is independent of DNS
entirely, so it survives the property being unset on the next machine, which is
the failure it was actually bought against — but it is not waiting on anyone
else's fix. As of 2026-09-11 the DNS name resolves and the alias runs alongside
it, which is the arrangement to keep: two independent ways to reach the
container, neither of which needs the other to be working.

**One thing left that is more workaround than it needs to be.**
`sshd-cmux.conf` sets `UsePAM no`, which is why the `node` account has to be
unlocked with `usermod -p '*'` — OpenSSH refuses a locked account for public key
auth when PAM is off. `UsePAM yes` is Debian's own default and would need
neither. It is left alone because the current arrangement is proven and swapping
an authentication path that works for one that is merely more standard is a poor
trade; but that is the reason, not an argument that this is better.

The VS Code extension still works — `customizations.vscode` is read when it
attaches — but nothing depends on it.

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
reads as losing somewhere to write code — and cmux is a terminal with no editor
in it. But `workspaceMount` is a bind: the files are on the Mac the whole time,
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

## cmux

cmux has no devcontainer support and does not need any. It has something better:
**`cmux ssh` is a first-class workspace type with a Linux-side daemon, and a
container running sshd is a remote host like any other.** cmux's own integration
suites attach to a Docker container this way — `tests_v2/test_ssh_remote_docker_forwarding.py`
and friends, with an ephemeral published port, a throwaway key and host key
checking off — so it is a tested configuration rather than a clever one.

That leaves two ways in, and they are not equivalent:

| | `up.sh` (`devcontainer exec`) | `cmux-attach.sh` (`cmux ssh`) |
|---|---|---|
| extra surface in the image | none | `openssh-server`, a published loopback port |
| `cmux` CLI inside the container | no | **yes** |
| terminal survives cmux restarting | no | yes, reconnects |
| sidebar metadata, sftp file drop | partial | yes |
| browser pane egress | the host's | the container's, so inside the firewall |
| `TERM` in the container | pinned to `xterm-256color` | the host's own, `COLORTERM` included |

Use `up.sh` for a shell. Use `cmux-attach.sh` for the workspace you actually
work in.

### How the ssh path works

`cmux ssh` probes the remote platform, uploads a release-pinned `cmuxd-remote`
binary verified against a SHA-256 manifest embedded in the app, and runs it over
stdio. **The daemon arrives over the SSH connection, not from the internet**,
which is why the firewall does not have to be opened to allow any of this.

It then reverse-forwards a TCP port — `ssh -N -R` — to an authenticated local
relay, installs a `cmux` wrapper at `~/.cmux/bin/cmuxd-remote`/`bin/cmux` on the
remote, prepends that to `PATH`, and pins `CMUX_SOCKET_PATH=127.0.0.1:<port>` in
the session. The relay port is per workspace. That is the whole trick: it is why
the `cmux` CLI works *from inside the container*.

What this repository adds for it:

- `openssh-server` in the image, and `sshd-cmux.conf` — cmux's own test fixture,
  minus root login and moved to port 2222. `AllowTcpForwarding yes` is
  load-bearing, not boilerplate: the reverse forward runs with
  `ExitOnForwardFailure=yes`, so refusing it fails the attach rather than
  degrading it.
- `start-sshd.sh`, reachable through `sudo` for the same reason
  `init-firewall.sh` is, and run from `postStartCommand` after it — the firewall
  flushes the tables, so anything holding a connection wants to start after it.
  Host keys are generated at first start rather than baked into the image, so two
  containers from one image do not share one.
- `-p 127.0.0.1::2222` in `runArgs`. No host port is named, so Docker picks a
  free one and several of these containers coexist; `cmux-attach.sh` finds it
  with `docker port`. Bound to loopback, so it is not on the network.
- `usermod --shell /bin/zsh node`. sshd reads the login shell out of
  `/etc/passwd` and ignores `ENV SHELL`, so without this an ssh session lands in
  bash — no `~/.zshrc`, so no `gh-auth.sh`, so no `gh`. Everything else names
  zsh explicitly and is unaffected.

Starting sshd unconditionally is safe: the port is loopback-only, password auth
is off, and no `authorized_keys` exists until `cmux-attach.sh` writes one. Until
you ask for it, it listens and refuses everything.

The key is a dedicated one, `~/.ssh/cmux-devcontainer`, not the one that talks to
GitHub — giving a loopback hop into a sandbox a key with any other reach is how a
sandbox stops being one. Agent forwarding is explicitly off for the same reason.
Host key checking is off because a host key is generated per container and the
port is a fresh ephemeral one each rebuild, so `known_hosts` could only ever
reject a container it had seen before on a port something else used. What bounds
this is the port being on loopback and the key being one the script made.

### The same sshd serves a second frontend

The Claude Code desktop app has an SSH environment — its documentation names dev
containers as a target — and it asks for exactly what `cmux-attach.sh` already
produces: a host, a port, and an identity file. `node@127.0.0.1`, whatever
`docker port` reports, and `~/.ssh/cmux-devcontainer`. It installs Claude Code
on the remote itself and uses the remote `/home/node/.claude`, which is the
volume. So the sshd here is not cmux-specific, and a native app is a second way
in rather than a different setup.

Two things to expect if you try it. There is no terminal panel in a remote
session, so it wants a terminal beside it rather than replacing one. And SSH
mode sets `CLAUDE_CODE_PROVIDER_MANAGED_BY_HOST=1` with its own
`ANTHROPIC_BASE_URL` and token, overriding provider configuration on the remote
— which costs nothing here, and would matter to a container pointed at Bedrock
or Vertex.

### What now works inside the container

- **`cmux notify`**, so a Claude Code hook can raise the ring and the sidebar
  badge directly rather than by printing an escape sequence. The relay is
  per-workspace, which is what lets cmux resolve which workspace a notification
  came from. Worth confirming once on the machine.
- `cmux workspace loading on`, `cmux read-screen`, `cmux send` — the CLI is
  relayed as a whole, not a subset.
- Terminals that survive cmux quitting, and reconnect on relaunch.

`install-cmux-hooks.sh`, run once per `~/.claude` volume from inside the
container, uses that to fix the one thing cmux gets wrong about a containerised
agent. It only draws its status pill for a process it recognises as an agent,
and it recognises the `claude` its own wrapper started on the host — not one in
here, which is a generic process to it. So the hooks drive the documented
workspace lane instead: `cmux workspace status set needs-attention` when Claude
asks something, `auto` when it stops. The sidebar row goes amber and back on its
own. Every hook is guarded on `command -v cmux`, so it is silent on the `up.sh`
path rather than an error every turn.

An OSC escape sequence still works too, and needs nothing installed:

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
belongs in `/home/node/.claude`, which is a persisted volume.

### What still does not work

**Native Claude session restore.** cmux's Claude Code integration is a wrapper
around the `claude` binary *on the host*, and its session records live in
`~/.cmuxterm/` there. A `claude` inside a container is not wrapped, so no
automatic `claude --resume <id>` on relaunch, no agent hibernation, and no AI
workspace naming. The manual equivalent, once per pane:

```bash
cmux surface resume set --shell '/path/to/repo/.devcontainer/cmux-attach.sh claude'
```

cmux keeps a socket-set command for manual restore until its prefix is approved
under Settings › Terminal › Resume Commands, which is deliberate on its part.

**Browser panes are inside the firewall.** cmux routes a remote workspace's
browser through a SOCKS5 proxy tunnelled over the daemon, so it egresses from the
container — and the allowlist blocks nearly everything. That is correct rather
than broken, but it means browsing happens in a local workspace, which cmux does
not force-proxy.

**`--transport mosh`.** Mosh needs inbound UDP in the 60000 range and the
firewall drops all UDP but DNS. Stay on the SSH transport.

### Other answers to the same problem

[cmux-devcontainer-bridge](https://github.com/zackey-heuristics/cmux-devcontainer-bridge)
is a Go daemon on the host listening on `127.0.0.1:8765`. A Claude Code hook in
the container `POST`s to `host.docker.internal:8765/notify`, and the bridge
execs `cmux notify` on the host. It exists for exactly the gap this file
describes, and it solves it from the other side: instead of getting a `cmux`
into the container, it gets the container's message out to the one on the host.

Not adopted here, for two reasons rather than one. Over `cmux ssh` the real CLI
is already in the container, authenticated per workspace, with nothing listening
on the host. And on the `up.sh` path the OSC hook above needs no daemon at all.
What the bridge adds over OSC is structured title/subtitle/body and the ability
to name a workspace by ID when the hook has no tty to write to — real, but
narrow.

Two things worth knowing before reaching for it anyway. `--token` is empty by
default, so anything in the sandbox that can reach the host gateway can drive
`cmux notify`; and nothing here would need opening for that, because
`init-firewall.sh` already accepts the host network in both directions. The exec
itself is safe by construction — `internal/notifier/cmux.go` hardcodes the
`notify` subcommand and passes values as argv rather than a shell string, and
the server uses a constant-time token compare and a body limit. It is one
release, one author and no stars, so build it from source; it has no third-party
dependencies, which makes that easy. Its example overlay assumes the cmux
devcontainer's router/sandbox compose split, which this container does not have.

What was worth taking from it outright is in `up.sh`:
`devcontainer exec --remote-env CMUX_WORKSPACE_ID=... CMUX_SURFACE_ID=...`.
Without it a hook inside the container has no idea which pane it belongs to.

[ccmanager](https://github.com/kbwo/ccmanager) is the other shape: it runs the
agent session inside the devcontainer as a first-class feature, with the manager
on the host. The trade is a session manager instead of a terminal — no browser
pane, no splits, no socket API, and none of the above.

### A different runtime under all of this

Docker is not the only way to get a Linux container on a Mac, and the one that
would suit this best is Apple's own — every container gets an address reachable
from the host, which deletes the published port and the `docker port` step from
`cmux-attach.sh` rather than adding to them. `NET_ADMIN` is supported, so the
firewall survives the move. What does not work yet is `build.dockerfile`, which
is what this repository uses. [APPLE.md](./APPLE.md) has the evidence, the
config to swap in, and the one upstream issue to watch. OrbStack is the
meanwhile option: it makes Docker faster without making it different, so nothing
here changes.

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

## Copying this into another project

Copy two directories and run two commands. Nothing in them names this project:

```bash
cp -R .devcontainer .cmux ../otherproject/
cd ../otherproject
./.devcontainer/build.sh
./.devcontainer/up.sh            # or the "Open sandbox" palette entry
```

Then once inside, per project, because both write to that project's own volume:

```bash
gh auth login
/workspace/.devcontainer/install-cmux-hooks.sh
```

**One folder name drives everything**, which is the only thing to get right.
`devcontainer.json` takes `name` and `image` from
`${localWorkspaceFolderBasename}`, `build.sh` derives the same tag with
`basename`, and `state.sh` keys the volume on the same string. So a checkout in
`~/src/pinpoint` is the container `pinpoint`, the image
`pinpoint-devcontainer:local`, the volume `claude-code-state-pinpoint` and the
hostname `pinpoint.adevcontainers.local`, with nothing written down anywhere.

Keep it lowercase and DNS-safe — it becomes a hostname and an image tag — and
**keep it distinct across projects**, because the volume name is that basename:
two checkouts both called `api` in different parents would share one Claude
config, history and `gh` login, silently. That is the one real hazard here.

The single optional edit is cosmetic: `.cmux/cmux.json` sets the palette
workspace's title. Left alone it says `proteolyzer`; the alternative is dropping
the key and letting cmux title the tab from the full path, which is uglier. Every
other value in that file is already derived — the commands run
`git rev-parse --show-toplevel`.

Two things this does not carry, both deliberate. The `agent`-label workflow wants
the GitHub App installed on the new repository, and the DNS name wants
`dns.domain` set once per *machine* rather than per project — see above, and both
are one-time rather than per-copy.

One edge worth knowing: run `build.sh` from a worktree rather than the checkout
root and the basename is the worktree's, so you get a separate image and volume.
That is consistent with `up.sh` from a worktree giving a separate container, and
it is usually what you want, but it is not what you want by accident.

### Converting an existing Docker devcontainer

A project already running under VS Code and Docker has the interesting part
already: its Claude config, session history, memory and shell history, sitting in
Docker's volume store. Bringing it over is `migrate` then `export`/`import`,
because the two runtimes keep separate stores and nothing is shared between them.

**First, read the old volume names off Docker rather than guessing**, and keep the
old `devcontainer.json` open long enough to see how it named its mounts:

```bash
docker volume ls | grep -i claude
```

The defaults in `state.sh` assume the folder basename —
`claude-code-config-pinpoint` and friends. Anthropic's own template keys them on
`${devcontainerId}` instead, which resolves to a hash, so a project that started
from that template has names like `claude-code-config-a1b2c3d4`. `migrate` looks
up by name and will otherwise report three misses and stop, which looks like "no
state to move" and is not.

`migrate` prints one line per volume — `will read …` or `no …; skipping` — and
reading those back is the check. A single skipped `gh` volume is the easy one to
miss: the fold succeeds, the container comes up, and the GitHub login is the one
thing that did not arrive.

Then, with the old container stopped — a volume attaches to one container at a
time, so this is a real requirement rather than tidiness:

```bash
cp -R ../proteolyzer/.devcontainer ../proteolyzer/.cmux .    # gives you state.sh

# fold the old three into one, inside Docker's store. Names as found above --
# set each of the three you actually have, gh included.
OLD_CONFIG=claude-code-config-a1b2c3d4 \
OLD_HISTORY=claude-code-bashhistory-a1b2c3d4 \
OLD_GH=claude-code-gh-a1b2c3d4 \
  RUNTIME=docker ./.devcontainer/state.sh migrate

# carry it across the runtime boundary
RUNTIME=docker    ./.devcontainer/state.sh export --with-credentials ~/s.tar.gz
RUNTIME=container ./.devcontainer/state.sh import --with-credentials ~/s.tar.gz

./.devcontainer/build.sh && REBUILD=1 ./.devcontainer/up.sh
rm ~/s.tar.gz
```

What arrives: the Claude login, the project's session history and memory, the
shell history, and — with `--with-credentials` — the GitHub token. Drop that flag
and everything else still comes, with `gh auth login` to run once on the far side.

Three things worth knowing. **The folder must keep its name**, because the new
volume is keyed on the basename and the container looks it up by that and nothing
else. **The old Docker volumes are left untouched**, so this is reversible by
deleting the new volume and reopening the project in VS Code. And **delete the
archive afterwards** — it holds the Claude credential whether or not you passed
`--with-credentials`, and the repository directory is often a synced folder.

### Skipping the build entirely

`build.sh` is the only per-project step that costs real time, and it does not
have to exist. `publish.sh` already pushes a multi-arch image to GHCR and takes
`IMAGE` for a neutral name:

```bash
IMAGE=ghcr.io/orhunkok/claude-devcontainer ./.devcontainer/publish.sh
```

Point `image:` at that tag instead of `${localWorkspaceFolderBasename}-devcontainer:local`
and a new project is `cp -R`, then `up.sh` — no build, no first-run wait. The
trade is direction of coupling: one image for every project means a Dockerfile
change is a publish plus a rebuild everywhere, rather than a local rebuild of the
one project you are working on. Worth it once the Dockerfile stops changing
weekly; not before.

## Moving to another machine

Most of it already travels, and it is worth being exact about which part does
not:

| | travels | how |
|---|---|---|
| the environment itself | yes | `Dockerfile` and `devcontainer.json`, in git |
| the workspace | yes | a bind mount, so wherever you cloned |
| the cmux commands | yes | `.cmux/cmux.json`, in git |
| the ssh key and host key | yes | `cmux-attach.sh` and `start-sshd.sh` make them if absent |
| **the three named volumes** | **no** | `state.sh` |

The volumes are the gap, and one of them matters: `/home/node/.claude` holds
settings, project state and the agent's memory directory. `state.sh export`
writes them to a tarball you copy across and `state.sh import` puts them back.
The `gh` volume is excluded unless you ask for it, because including it writes a
GitHub token into that tarball in plaintext and `gh auth login` is one command
on the far side. Volume names are keyed on the directory basename, so the
checkout has to be named the same over there — the same constraint
`devcontainer.json` already documents.

**The macOS-only parts are the frontend, not the environment.** cmux does not
run on Linux and neither does Apple `container` — but neither of them is the
environment. The image is OCI, so Docker, podman, containerd and Apple
`container` all take it and none of them cares which one built it; the runtime
underneath is interchangeable and picking a Mac-native one costs nothing here.
What is runtime-specific is the two attach scripts, and that is why there are
two: `cmux-attach.sh` is the nice thing on a Mac, `up.sh` is the one that still
works on a Linux box, and neither is the source of truth.

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
