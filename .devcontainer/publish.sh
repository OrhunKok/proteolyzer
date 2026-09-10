#!/usr/bin/env bash
# Publish this container's image to GHCR, multi-arch, so another machine pulls
# it instead of building it.
#
# The image is the portable unit -- OCI, so Docker, podman, containerd and Apple
# `container` all take it, and none of them cares which one built it. What makes
# a machine change slow is rebuilding from the Dockerfile, not the format. This
# removes that.
#
# Runs on the host. Needs `gh` logged in with `write:packages`.
#
#   ./.devcontainer/publish.sh              build both arches and push
#   ./.devcontainer/publish.sh --dry-run    print the tags and stop
#
# Afterwards, `image:` in devcontainer.json can point at the tag it prints,
# instead of `build:`. Keep `features` where it is: both the devcontainer CLI and
# adevcontainer derive an image from base + features at create time, so the
# Python feature is unaffected by the switch.
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
dry_run=0
[ "${1:-}" = "--dry-run" ] && dry_run=1

slug="$(git -C "$repo" remote get-url origin 2>/dev/null \
    | sed -E 's#(git@github\.com:|https://github\.com/)##; s#\.git$##')"
if [ -z "$slug" ]; then
    echo "publish.sh: no github origin to derive the image name from." >&2
    exit 1
fi

# GHCR wants the owner lowercased; the repository part may keep its case but
# there is no reason to find out the hard way.
#
# The default names the image after this repository, which is wrong the moment a
# second project wants it -- and nothing in the Dockerfile is specific to this
# one. Set IMAGE to a neutral name and every project can pin the same image:
#
#   IMAGE=ghcr.io/orhunkok/claude-devcontainer ./.devcontainer/publish.sh
image="${IMAGE:-ghcr.io/$(echo "$slug" | tr '[:upper:]' '[:lower:]')-devcontainer}"

# Tag on the content of .devcontainer rather than the git SHA: the image only
# changes when the thing that builds it does, so an unrelated commit does not
# invalidate a pull, and a config can pin a tag that means something. `ls-files`
# reads the index, so an edit you have not staged does not move the tag -- stage
# or commit before publishing, or you will push new content under an old name.
digest="$(git -C "$repo" ls-files -s .devcontainer | git hash-object --stdin | cut -c1-12)"

echo "publish.sh: $image:$digest"
echo "publish.sh: $image:latest"
[ "$dry_run" -eq 1 ] && exit 0

for tool in docker gh; do
    command -v "$tool" >/dev/null 2>&1 || { echo "publish.sh: $tool is not on PATH." >&2; exit 1; }
done

# buildx cross-builds the arch this machine is not, which is the whole point:
# one push serves an Apple silicon laptop and an x86 box.
if ! docker buildx version >/dev/null 2>&1; then
    echo "publish.sh: docker buildx is required for a multi-arch build." >&2
    exit 1
fi

gh auth token | docker login ghcr.io -u "$(gh api user --jq .login)" --password-stdin

docker buildx build \
    --platform linux/amd64,linux/arm64 \
    --tag "$image:$digest" \
    --tag "$image:latest" \
    --push \
    "$repo/.devcontainer"

echo
echo "publish.sh: pushed. GHCR makes a new package private by default;"
echo "publish.sh: a machine that pulls it needs a token with read:packages."
