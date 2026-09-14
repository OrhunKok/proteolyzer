#!/usr/bin/env bash
# Build the image `devcontainer.json` names.
#
# This exists because the config is `image:` rather than `build:` — which is what
# lets Apple `container` run it at all, since adevcontainer rejects Dockerfile
# builds. The trade is that the Dockerfile stops being built implicitly: change
# it and run this before the change takes effect.
#
# Runs on the host, on Apple `container`. What it produces is ordinary OCI, so
# any runtime would take it -- `container build` is simply the one installed here.
#
#   ./.devcontainer/build.sh
#   IMAGE=my-devcontainer:local ./.devcontainer/build.sh
set -euo pipefail

# These drive macOS-side tooling -- `container`, `adevcontainer`, Docker on the
# Mac -- so running one *inside* the container is a mistake worth naming. Left
# uncaught the symptom is a tool "not on PATH" plus an invitation to
# `brew install` it, on Linux, which sends you somewhere with no exit.
if [ "$(uname -s)" = Linux ]; then
    printf '%s\n' \
        "${0##*/}: this runs on the Mac, not inside the container." \
        "${0##*/}: \`exit\` back to the host first, or open a terminal there." >&2
    exit 1
fi

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Read out of devcontainer.json rather than derived independently, so the image
# this builds and the image adevcontainer starts cannot drift apart. The one
# substitution that file uses is applied here too. IMAGE still overrides.
config="$repo/.devcontainer/devcontainer.json"
declared="$(sed -n 's/^[[:space:]]*"image"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' "$config" | head -1)"
tag="${IMAGE:-${declared//\$\{localWorkspaceFolderBasename\}/$(basename "$repo")}}"

if [ -z "$tag" ]; then
    echo "build: no \"image\" in $config." >&2
    exit 1
fi

# An OCI reference has to be lowercase, and a folder name does not. streamlit-DO-MS
# derives streamlit-DO-MS-devcontainer:local, which Apple `container` rejects as
# `invalid reference` -- an error that names neither the cause nor the fix, and
# does not say which of the two files to change. Both read the same declaration
# now, so there is only one place to change, and this says so before the build.
case "$tag" in
    *[A-Z]*)
        lower="$(printf '%s' "$tag" | tr '[:upper:]' '[:lower:]')"
        echo "build: \"$tag\" is not a legal image reference -- they must be lowercase," >&2
        echo "build: and this folder's name is not. Set it explicitly in" >&2
        echo "build: $config:" >&2
        echo "build:" >&2
        echo "build:     \"image\": \"$lower\"," >&2
        echo "build:" >&2
        echo "build: then rerun. Both this script and adevcontainer read that line," >&2
        echo "build: so they cannot disagree about which image is meant." >&2
        exit 1
        ;;
esac

if ! command -v container >/dev/null 2>&1; then
    echo "build: Apple \`container\` is not on PATH." >&2
    echo "build:   brew install --cask container      (macOS 26+)" >&2
    exit 1
fi

# Apple `container` rejects a build request over roughly 16 KiB
# (apple/container#735). Well clear of it today; worth saying before someone
# doubles the Dockerfile and gets `Stream unexpectedly closed` instead.
size="$(wc -c < "$repo/.devcontainer/Dockerfile")"
if [ "$size" -gt 14336 ]; then
    echo "build: warning: Dockerfile is ${size} bytes; Apple container's build" >&2
    echo "build: request limit is around 16 KiB and failures there are opaque." >&2
fi

# The identity this Mac commits with, baked into the image as the default for
# every project built from it -- otherwise `git config --global user.email` is a
# step in each one, and the symptom of forgetting arrives much later as
# `Author identity unknown` on a commit. Read from the full chain rather than
# `--global`, so a repository-local identity counts: the question is what you
# actually commit as, not where you wrote it down.
git_name="$(git -C "$repo" config user.name 2>/dev/null || true)"
git_email="$(git -C "$repo" config user.email 2>/dev/null || true)"
if [ -z "$git_email" ]; then
    echo "build: this Mac has no git identity, so the image will carry no" >&2
    echo "build: default one. Set \`git config --global user.email\` and rebuild," >&2
    echo "build: or set it inside each container." >&2
fi

echo "build: container build -t $tag"
container build \
    --build-arg TZ="${TZ:-America/New_York}" \
    --build-arg CLAUDE_CODE_VERSION=latest \
    --build-arg GIT_DELTA_VERSION=0.18.2 \
    --build-arg ZSH_IN_DOCKER_VERSION=1.2.0 \
    --build-arg GIT_USER_NAME="$git_name" \
    --build-arg GIT_USER_EMAIL="$git_email" \
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
stale="$(container image list 2>/dev/null | awk '/adev-/ {print $1":"$2}' || true)"
if [ -n "$stale" ]; then
    echo "build: dropping derived images so the new base is actually used:"
    printf '%s\n' "$stale" | sed 's/^/build:   /'
    printf '%s\n' "$stale" | while read -r image; do
        container image delete "$image" >/dev/null 2>&1 || \
            echo "build: could not delete $image; delete it by hand" >&2
    done
else
    # Not fatal, but said out loud: the column layout of `container image list`
    # is not a contract, and a silently missed purge presents as a change that
    # did not apply -- this exact bug.
    echo "build: no derived adev- images found. If a change does not take" >&2
    echo "build: effect, check \`container image list\` for one and delete it." >&2
fi

echo "build: next: REBUILD=1 ./.devcontainer/up.sh"
