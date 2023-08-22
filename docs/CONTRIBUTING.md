# Contributing to CHAMPAGNE

TODO -- describe gitflow, require PRs...

### Use pre-commit hooks

Pre-commit can automatically format your code, check for spelling errors, etc. every time you commit.

Install [pre-commit](https://pre-commit.com/#installation) if you haven't already,
then run `pre-commit install` to install the hooks specified in `.pre-commit-config.yaml`.
Pre-commit will run the hooks every time you commit.

### Versions

Increment the version number following semantic versioning[^3] in `src/champagne/VERSION`.

[^3]: semantic versioning guidelines https://semver.org/

### Changelog

Keep the changelog up to date with any user-facing changes in `docs/CHANGELOG.md`.