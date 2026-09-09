# The sandbox, without VS Code

Two reasons, and the order matters when weighing a replacement. The first is
that the environment is a file in git, so another machine is a clone rather than
an afternoon. The second is that a firewall stands between Claude Code and
everything that is not GitHub, PyPI, npm or the Anthropic API — see
`init-firewall.sh`, which is why `runArgs` asks for `NET_ADMIN`.

Most things that claim to replace this replace the second reason only. See
[Moving to another machine](#moving-to-another-machine).

The container is defined by `devcontainer.json` and nothing about it needs an
editor:

```bash
npm install -g @devcontainers/cli    # once
./.devcontainer/up.sh                # a zsh in the container
./.devcontainer/up.sh claude         # straight into Claude Code
./.devcontainer/cmux-attach.sh       # the same container as a cmux workspace
REBUILD=1 ./.devcontainer/up.sh      # discard the container and build again
```

`up.sh` is `devcontainer up` followed by `devcontainer exec`, which is the whole
of what the Dev Containers extension was doing to the container itself. The
extension still works — `customizations.vscode` is read when VS Code attaches —
but nothing depends on it.

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

`~/.config/gh` is a named volume keyed on the directory name, so that login
survives a rebuild the way the command history and the Claude config do.
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

`streamlit-DO-MS` and `decoder` carry their own copy of this container. Nothing
here reaches into them: the change is a handful of files and a `gh auth login`,
applied in each repository by that repository, which is the same rule as
everything else on the account.
