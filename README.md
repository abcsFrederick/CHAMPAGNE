# CHAMPAGNE 🍾

**CH**rom**A**tin i**M**muno **P**recipit**A**tion sequencin**G** a**N**alysis pip**E**line

🚧 **This project is under active development. It is not yet ready for production use.** 🚧

## Getting started

TODO

## Usage

### champagne CLI

You can run champagne from the command line.
The CLI includes helper steps for execution on supported
high performance computing clusters including Biowulf and FRCE.

Install the champagne CLI:

```sh
cd CHAMPAGNE
pip3 install .
```

Run the test dataset using the test profile:

```sh
champagne run -profile test,singularity
```

or explicitly specify the output directory and input:

```sh
champagne run -profile singularity --outdir results/test --input assets/samplesheet_test.csv
```

Launch a stub run to view the steps that will run without performing the full analysis.

```sh
champagne run -profile test -stub
```

### nextflow pipeline

You can run the nextflow pipeline directly by specifying this GitHub repo.
You will need nextflow and either singularity or docker installed.

```sh
nextflow run CCBR/CHAMPAGNE -profile test,singularity
```

You can specify a specific version, tag, or branch with `-r`:

```sh
nextflow run CCBR/CHAMPAGNE -r v1.0.0 -profile test,singularity
```

## Help & Contributing

Come across a **bug**? Open an [issue](https://github.com/CCBR/CHAMPAGNE/issues) and include a minimal reproducible example.

Have a **question**? Ask it in [discussions](https://github.com/CCBR/CHAMPAGNE/discussions).

Want to **contribute** to this project? Check out the [contributing guidelines](docs/CONTRIBUTING.md).

## References

This repo was originally generated from the
[CCBR Nextflow Template](https://github.com/CCBR/CCBR_NextflowTemplate).
The template takes inspiration from nektool[^1] and the nf-core template.
If you plan to contribute your pipeline to nf-core, don't use this template --
instead follow nf-core's instructions[^2].

[^1]: nektool https://github.com/beardymcjohnface/nektool
[^2]: instructions for nf-core pipelines https://nf-co.re/docs/contributing/tutorials/creating_with_nf_core
