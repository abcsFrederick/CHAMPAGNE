# CHAMPAGNE 🍾

**CH**rom**A**tin i**M**muno **P**recipit**A**tion sequencin**G** a**N**alysis pip**E**line

 🚧**This project is under active development. It is not yet ready for production use.** 🚧

## Getting started

TODO

## Usage

Install the tool in edit mode:

```sh
pip3 install -e .
```

Run the example

```sh
champagne run --input "Hello world"
```

## Developer notes

### Use pre-commit hooks

Pre-commit can automatically format your code, check for spelling errors, etc. every time you commit.

Install [pre-commit](https://pre-commit.com/#installation) if you haven't already,
then run `pre-commit install` to install the hooks specified in `.pre-commit-config.yaml`.
Pre-commit will run the hooks every time you commit.

### Versions

Increment the version number following semantic versioning[^3] in `src/CHAMPAGNE/VERSION`

[^3]: semantic versioning guidelines https://semver.org/

### Changelog

Keep the changelog up to date with any user-facing changes in `docs/CHANGELOG.md`

## References

This repo was originally generated from the [CCBR Nextflow Template](https://github.com/CCBR/CCBR_NextflowTemplate).
The template takes inspiration from nektool[^1] and the nf-core template. If you plan to contribute your pipeline to nf-core, don't use this template -- instead follow nf-core's instructions[^2].

[^1]: nektool https://github.com/beardymcjohnface/nektool
[^2]: instructions for nf-core pipelines https://nf-co.re/docs/contributing/tutorials/creating_with_nf_core
