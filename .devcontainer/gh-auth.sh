# Make `gh` work inside this devcontainer.
#
# The container inherits a GitHub credential from the host through the VS Code
# credential helper, which is why `git push` works. `gh` does not read git
# credential helpers, so it reports itself as logged out with the credential
# sitting right there.
#
# The token is fetched per invocation rather than exported, so it is not sitting
# in the environment of every child process.
#
# Sourced from ~/.zshrc by postCreateCommand. This is the whole of what used to
# be a shared repository, a mounted checkout, an installer and two session hooks
# -- see DECISIONS.md.

gh() {
    local token
    token="$(printf 'protocol=https\nhost=github.com\n\n' \
        | git credential fill 2>/dev/null \
        | sed -n 's/^password=//p')"

    if [ -z "$token" ]; then
        echo "gh-auth: no GitHub credential available from git credential fill." >&2
        echo "gh-auth: is the host's credential helper forwarded to this container?" >&2
        return 1
    fi

    GH_TOKEN="$token" command gh "$@"
}
