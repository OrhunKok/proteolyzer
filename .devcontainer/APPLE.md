# Apple's stack, and where this container stops

Two questions, both about the same boundary. First: can this sandbox run on
Apple's own container runtime instead of Docker. Second: can Apple's on-device
AI stack, Core AI, be used from inside it. Nothing below was run — this machine
has neither macOS nor any of the tools. It is read off the projects' own
specifications and PyPI metadata, cited so the claims can be rechecked when
they move.

## Running the sandbox on Apple `container`

**Verdict: the runtime is ready, the devcontainer layer is one feature short.**

The appeal is real and it is not just tidiness. Under Docker this container is
reached through a published loopback port that Docker picks at random, which is
why `cmux-attach.sh` has a `docker port` step and a
`-p 127.0.0.1::2222` in `runArgs`. Apple's networking doc removes all of that:

> Every container gets an IP address on its network, always reachable by that IP
> from the host and from other containers on the same network (find it with
> `container inspect <name>`).

And with `domain = "test"` under `[dns]` in `~/.config/container/config.toml`,
every container registers as `<name>.test` and macOS can be pointed at the same
resolver. So the attach is `cmux ssh node@proteolyzer.test` — no published port,
no port discovery, no loopback binding, and `cmux-attach.sh` gets shorter rather
than longer. That is a better answer than OrbStack, which makes Docker faster
without making it different.

## What already works

`--cap-add NET_ADMIN` is supported, which is the one that decides it —
`init-firewall.sh` is the reason this container exists, and it needs
`NET_ADMIN` for `iptables` and `ipset`. Apple's `docs/runtime-configuration.md`
gives the literal invocation:

```bash
container run --cap-add NET_ADMIN alpine ip link set lo down
```

`CAP_NET_RAW` is already in the default capability set, so the second half of
this repository's `runArgs` is redundant there rather than blocked.

[`wcgomes/apple-devcontainers`](https://github.com/wcgomes/apple-devcontainers)
(`adevcontainer`) is a Swift CLI that reads `devcontainer.json` and drives Apple
`container` directly — not a Docker shim. Its `specs/core.md` admits `capAdd` as
a first-class property on a "typed capability path", plus `mounts` (bind and
volume), `containerEnv`, `remoteUser`, `workspaceFolder`, the full lifecycle
hook set with `waitFor`, and `customizations.vscode`. It has a real Features
runner, so `ghcr.io/devcontainers/features/python:1` installs.

## What blocks it

**`build.dockerfile` hard-errors.** `specs/core.md` lists "unrepresentable
configuration-source selectors (`build`/legacy Dockerfile and Compose until
separately supported)" among the hard rejections. This repository's config is
`build.dockerfile`, so `adevcontainer up` refuses it outright. The capability is
there and only the config surface is missing: the Features runner already
generates a Dockerfile and shells out to `container build` to derive its image.

**`-p` in `runArgs` hard-errors.** ADR 0003 fails closed on "first-class
smuggling via runArgs (`-e`, `-u`, `-w`, `-p`, `-v`, …)". Which costs nothing,
because the container has its own address and does not need publishing.
`--cap-add` is not in that list and is passed through.

**`workspaceMount` is unsupported.** There is an implicit bind of the host
workspace root to the container workspace folder instead, so the line comes out.

**Config discovery is fixed** to `.devcontainer/devcontainer.json` then
`.devcontainer.json`, with no way to name a different file. So an Apple variant
cannot sit beside the Docker one — it is a swap, not an addition. That is the
reason this file is prose and not a second `devcontainer.json`.

## A runtime papercut worth knowing before it bites

Apple `container` does not seed a named volume from the image path it covers —
[apple/container#729](https://github.com/apple/container/issues/729), open since
October 2025, with a reproduction in the issue body. The `Dockerfile` here
depends on exactly the opposite, and says so:

> a named volume is initialised from the image path it covers, ownership
> included, so a path the image does not have comes up root-owned and unwritable
> by node

Under Apple `container` all three volumes — `/commandhistory`,
`/home/node/.claude`, `/home/node/.config/gh` — come up empty and root-owned,
and `sudo` here is restricted to two named scripts. That is the decoder
`Permission denied` again, from the other direction. The fix is a third script
on the same sudo arrangement, chowning those three paths, ahead of the firewall
in `postStartCommand`.

## The config, for when `build` lands

Swap `.devcontainer/devcontainer.json` for this, having built the image first:

```bash
container build -t proteolyzer-dev:latest .devcontainer
```

```jsonc
{
  "name": "Claude Code Sandbox",
  "image": "proteolyzer-dev:latest",

  // Was runArgs. NET_RAW is in Apple's default capability set, and -p is both
  // rejected and unnecessary -- the container has its own IP.
  "capAdd": ["NET_ADMIN"],

  "features": {
    "ghcr.io/devcontainers/features/python:1": {
      "version": "3.14",
      "installTools": true
    }
  },
  "remoteUser": "node",
  "mounts": [
    "source=claude-code-bashhistory-${localWorkspaceFolderBasename},target=/commandhistory,type=volume",
    "source=claude-code-config-${localWorkspaceFolderBasename},target=/home/node/.claude,type=volume",
    "source=claude-code-gh-${localWorkspaceFolderBasename},target=/home/node/.config/gh,type=volume"
  ],
  "containerEnv": {
    "CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC": "true",
    "DISABLE_TELEMETRY": "true",
    "NODE_OPTIONS": "--max-old-space-size=4096",
    "CLAUDE_CONFIG_DIR": "/home/node/.claude",
    "POWERLEVEL9K_DISABLE_GITSTATUS": "true"
  },

  // workspaceMount is gone: the implicit bind covers it.
  "workspaceFolder": "/workspace",

  "postCreateCommand": "grep -q gh-auth.sh ~/.zshrc || echo 'source /workspace/.devcontainer/gh-auth.sh' >> ~/.zshrc",
  // fix-volume-perms.sh first, because of apple/container#729 above.
  "postStartCommand": "sudo /usr/local/bin/fix-volume-perms.sh && sudo /usr/local/bin/init-firewall.sh && sudo /usr/local/bin/start-sshd.sh",
  "waitFor": "postStartCommand"
}
```

The attach loses its whole middle section — no `docker port`, no published port,
no loopback:

```bash
cmux ssh "node@proteolyzer.test" \
    --identity ~/.ssh/cmux-devcontainer \
    --name proteolyzer \
    --no-forward-agent \
    --command "cd /workspace"
```

## How much to trust it

`apple/container` is 49.8k stars and pushed daily; the runtime is not the risk.
`adevcontainer` is 11 stars, one fork and one author, five weeks old at the time
of writing — but unusually disciplined for that age: ADRs, requirement specs
with scenarios, a compatibility-degradation report with stable codes, and an
opt-in `ADEVCONTAINER_STRICT_COMPATIBILITY=1`. Young, not sloppy.

So: stay on Docker, and watch one thing —
[apple/container#2112](https://github.com/apple/container/issues/2112), where
`adevcontainer`'s author offers it as a `container` plugin. When `build` support
lands, this is a swap of one file and one build command.

## Core AI, and why it stays on the host

[Core AI](https://developer.apple.com/documentation/coreai) is Apple's on-device
stack after Core ML: its own IR, models exported as standalone `.aimodel` files,
and a compiler and runtime behind them.
[`apple/coreai-torch`](https://github.com/apple/coreai-torch) lowers a
`torch.export.ExportedProgram` into that IR;
[`apple/coreai-models`](https://github.com/apple/coreai-models) carries export
recipes, a Swift runtime package, and agent skills.

**Running a model is macOS/iOS 27 and Xcode 27**, so that half is host-only and
not arguable — the same boundary as everything else in this file.

**Converting one could in principle happen in here, and cannot today**, for two
reasons worth writing down because both are silent until you try:

- `coreai-torch` is `py3-none-any`, but it hard-depends on `coreai-core`, which
  declares `requires_python = ">=3.10,<3.14"`. `devcontainer.json` pins the
  Python feature to **3.14**. Out of range on the version this container is
  built around.
- `coreai-core` has only ever published `macosx_26_0_arm64`, `manylinux1_x86_64`
  and `manylinux_2_34_x86_64` wheels — across both releases, **no Linux
  aarch64**. On an Apple silicon Mac this container is arm64 Linux, so there is
  no wheel for it at any Python version. That one is not fixable from this side.

The x86_64 manylinux wheels say Apple means the export path to be usable off a
Mac, so ARM Linux is plausibly a matter of time. Until then: on the host.

The agent skills in `coreai-models` (`working-with-coreai`, `model-authoring`,
`model-compression-exploration`) follow the runtime for the same reason. They
drive `coreai-torch` and `coreai-opt` against a Mac, so they belong to a Claude
Code running on the host, not to the one in here — even though
`/home/node/.claude` is a volume and would happily hold them.
