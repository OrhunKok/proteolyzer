# Make `gh` and `git` work inside this devcontainer, however it was started.
#
# There are two credential paths and only one of them is ever available, decided
# by what started the container:
#
#   VS Code forwards a credential from the host through its own git credential
#   helper, which is an RPC pipe back to the editor process. That is why `git
#   push` used to work for free. `gh` does not read git credential helpers, so it
#   reported itself as logged out with the credential sitting right there --
#   which is what the wrapper below is for.
#
#   `devcontainer up` forwards nothing: no editor, no pipe, so `git credential
#   fill` comes back empty and both tools are dead. The container has to hold its
#   own credential instead. `gh auth login` once, into the /home/node/.config/gh
#   volume so it survives a rebuild, and `gh auth setup-git` lends the same
#   credential to git -- one login, both tools, no host involved.
#
# Cheapest test first, and by file rather than by `gh auth status`, which spends a
# network round trip on every shell start.
#
# Sourced from ~/.zshrc by postCreateCommand. This is the whole of what used to
# be a shared repository, a mounted checkout, an installer and two session hooks
# -- see DECISIONS.md.

if [ -n "${GH_TOKEN:-${GITHUB_TOKEN:-}}" ]; then
    : # A token is already in the environment; gh reads it itself.

elif [ -s "${GH_CONFIG_DIR:-$HOME/.config/gh}/hosts.yml" ]; then
    # gh holds its own credential. Lend it to git as well, once per rebuild, and
    # only when nothing else has claimed the helper -- a second helper beside a
    # forwarded one is a coin toss over which answers.
    if ! git config --global --get-regexp '^credential\.' >/dev/null 2>&1; then
        gh auth setup-git 2>/dev/null || true
    fi

else
    # Nothing of our own: borrow the host's, fetched per invocation rather than
    # exported, so it is not in the environment of every child process.
    gh() {
        # `auth` has to reach the real gh, or `gh auth login` -- the fix this
        # very function tells you to run -- would be swallowed by the failure.
        if [ "${1:-}" = "auth" ]; then
            command gh "$@"
            return
        fi

        local token
        token="$(printf 'protocol=https\nhost=github.com\n\n' \
            | git credential fill 2>/dev/null \
            | sed -n 's/^password=//p')"

        if [ -z "$token" ]; then
            echo "gh-auth: no credential in ~/.config/gh, and none forwarded to" >&2
            echo "gh-auth: git credential fill either. Run: gh auth login" >&2
            return 1
        fi

        GH_TOKEN="$token" command gh "$@"
    }
fi
