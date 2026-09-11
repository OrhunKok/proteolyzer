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

# Derived from the folder name, not written down, so this script is the same in
# every project. devcontainer.json builds its `image` from
# ${localWorkspaceFolderBasename} and lands on the same string, and state.sh
# keys its volume on the same basename -- one source of truth, three files, no
# edit when you copy the directory. IMAGE still overrides for a one-off build.
tag="${IMAGE:-$(basename "$repo")-devcontainer:local}"

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

echo "build: $tag is ready."

# adevcontainer does not run this image directly. Because `features` is set, it
# derives one -- base plus features -- and reuses it under a tag hashed from the
# *config*, which does not include the base image. So rebuilding the base changes
# nothing: the container keeps starting from a derived image built before your
# edit, and the edit appears to have had no effect at all. A firewall change cost
# an hour to that, and the log calls it "Reusing features image adev-...", which
# reads like progress.
#
# Dropping the derived images makes it re-derive from what was just built. The
# cost is re-fetching the feature, which is seconds.
if [ "$runtime" = container ]; then
    stale="$(container image list 2>/dev/null | awk '/adev-/ {print $1":"$2}' || true)"
    if [ -n "$stale" ]; then
        echo "build: dropping derived images so the new base is actually used:"
        printf '%s\n' "$stale" | sed 's/^/build:   /'
        printf '%s\n' "$stale" | while read -r image; do
            container image delete "$image" >/dev/null 2>&1 || \
                echo "build: could not delete $image; delete it by hand" >&2
        done
    else
        # Not fatal, but said out loud: the column layout of
        # `container image list` is not a contract, and a silently missed purge
        # presents as a change that did not apply -- this exact bug.
        echo "build: no derived adev- images found. If a change does not take" >&2
        echo "build: effect, check \`container image list\` for one and delete it." >&2
    fi
fi

echo "build: next: REBUILD=1 ./.devcontainer/up.sh"
