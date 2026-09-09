# The sandbox, without VS Code

Claude Code runs in this container so that a firewall stands between it and
everything that is not GitHub, PyPI, npm or the Anthropic API — see
`init-firewall.sh`, which is the reason `runArgs` asks for `NET_ADMIN`.

The container is defined by `devcontainer.json` and nothing about it needs an
editor. Start it from a terminal:

```bash
npm install -g @devcontainers/cli   # once
./.devcontainer/up.sh               # a zsh in the container
./.devcontainer/up.sh claude        # straight into Claude Code
REBUILD=1 ./.devcontainer/up.sh     # discard the container and build again
```

`up.sh` is `devcontainer up` followed by `devcontainer exec`, which is the whole
of what the Dev Containers extension was doing to the container itself. The
extension is still supported — `customizations.vscode` is read when VS Code
attaches — but nothing depends on it any more.

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
old forwarding path is still there as a fallback for a VS Code session.

One edge: `gh auth setup-git` is skipped when a credential helper is already
configured, because two helpers is a coin toss over which answers. If you attach
VS Code to a container and then go back to the terminal, its helper is left
behind pointing at a pipe that is gone. `git config --global --unset-all
credential.helper` clears it, or rebuild.

**Formatting on save.** `editor.formatOnSave` and the eslint fixer in
`customizations.vscode` do nothing without the editor. `make lint` and
`.pre-commit-config.yaml` are what actually enforce formatting here, and they run
in the container either way.

Three entries in `init-firewall.sh` — `marketplace.visualstudio.com`,
`vscode.blob.core.windows.net`, `update.code.visualstudio.com` — exist so the VS
Code server and its extensions can install themselves inside the container. They
are dead weight once nothing attaches, and dropping them narrows the allowlist by
three. They are left in because they cost nothing and are the difference between
VS Code working and VS Code hanging on the day you want it back.

## cmux

cmux is a terminal, not a container manager. It replaces the editor half of what
VS Code was doing and none of the container half — there is no devcontainer
support in it, and there is no point looking for it: `devcontainer.json` stays,
`up.sh` drives it, and cmux is what the shell is displayed in.

`.cmux/cmux.json` in the repository root is read by cmux when a workspace is
opened here, and it registers the three commands above in the command palette
(⌘⇧P). Nothing else uses that file.

**Session restore.** cmux restores panes and working directories on relaunch, and
for an agent it runs the agent's native resume command — but that integration is
a wrapper around the `claude` binary *on the host*, and this `claude` is inside a
container where cmux cannot see it. What cmux can restore is the way back in:

```bash
cmux surface resume set --shell '/path/to/repo/.devcontainer/up.sh claude'
```

Run in the pane you want it attached to. cmux keeps a socket-set command for
manual restore until the prefix is approved under Settings › Terminal › Resume
Commands, which is deliberate on its part and worth doing once.

Losing the host wrapper also means no agent hibernation and no AI workspace
naming for a containerised session. Both are conveniences; the firewall is not.

**Notifications.** cmux raises its ring and sidebar badge off OSC 9/99/777, and
those pass through `devcontainer exec` to the terminal like any other escape
sequence, so a bell from inside the container still lands. `cmux notify` does
not — it is a host binary. To drive notifications from Claude Code's hooks,
write the sequence rather than calling the CLI, in the settings.json inside
`/home/node/.claude` (which is a persisted volume):

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

`> /dev/tty` is the part that matters: hook stdout is captured by Claude Code, so
a sequence merely printed never reaches the terminal. Confirm it on the machine
before relying on it.

**If what you wanted was devcontainers managed for you**, cmux is the wrong
shape and [ccmanager](https://github.com/kbwo/ccmanager) is the right one: it
keeps the session manager on the host and runs the agent session inside the
devcontainer as a first-class feature. The trade is a session manager instead of
a terminal — no browser pane, no splits, no socket API.

## Other repositories

`streamlit-DO-MS` and `decoder` carry their own copy of this container. Nothing
here reaches into them: the change is four files and a `gh auth login`, applied
in each repository by that repository, which is the same rule as everything else
on the account.
