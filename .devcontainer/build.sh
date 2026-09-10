#!/usr/bin/env bash
# Build the image `devcontainer.json` names.
#
# This exists because the config is `image:` rather than `build:` — which is what
# lets Apple `container` run it at all, since adevcontainer rejects Dockerfile
# builds. The trade is that the Dockerfile stops being built implicitly: change
# it and run this before the change takes effect.
#
# Runs on the host. Uses Apple `container` when it is installed, Docker
# otherwise; the image is ordinary OCI either way, so the one built by either
# runs under both.
#
#   ./.devcontainer/build.sh
#   IMAGE=my-devcontainer:local ./.devcontainer/build.sh
set -euo pipefail

# These drive macOS-side tooling -- cmux, `container`, `adevcontainer`, Docker
# on the Mac -- so running one *inside* the container is a mistake worth naming.
# Left uncaught the symptom is "cmux is not on PATH" plus an invitation to
# `brew install` it, on Linux, which sends you somewhere with no exit.
if [ "$(uname -s)" = Linux ]; then
    printf '%s\n' \
        "${0##*/}: this runs on the Mac, not inside the container." \
        "${0##*/}: \`exit\` back to the host first, or use another cmux tab." >&2
    exit 1
fi

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tag="${IMAGE:-proteolyzer-devcontainer:local}"

if command -v container >/dev/null 2>&1; then
    runtime=container
elif command -v docker >/dev/null 2>&1; then
    runtime=docker
else
    echo "build: no container runtime found." >&2
    echo "build:   brew install --cask container      (Apple, macOS 26+)" >&2
    echo "build:   or start Docker Desktop / OrbStack" >&2
    exit 1
fi

# Apple `container` rejects a build request over roughly 16 KiB
# (apple/container#735). Well clear of it today; worth saying before someone
# doubles the Dockerfile and gets `Stream unexpectedly closed` instead.
size="$(wc -c < "$repo/.devcontainer/Dockerfile")"
if [ "$runtime" = container ] && [ "$size" -gt 14336 ]; then
    echo "build: warning: Dockerfile is ${size} bytes; Apple container's build" >&2
    echo "build: request limit is around 16 KiB and failures there are opaque." >&2
fi

echo "build: $runtime build -t $tag"
"$runtime" build \
    --build-arg TZ="${TZ:-America/New_York}" \
    --build-arg CLAUDE_CODE_VERSION=latest \
    --build-arg GIT_DELTA_VERSION=0.18.2 \
    --build-arg ZSH_IN_DOCKER_VERSION=1.2.0 \
    -t "$tag" \
    "$repo/.devcontainer"

echo "build: $tag is ready. Next: ./.devcontainer/up.sh"
