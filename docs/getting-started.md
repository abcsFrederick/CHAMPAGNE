# Getting Started with CHAMPAGNE

TODO intro paragraph

## Installation

CHAMPAGNE is installed on the Biowulf and FRCE clusters as part of the
`ccbrpipeliner` module.
If you'd like to run the pipeline in a different execution environment,
take a look at [how to run the nextflow pipeline directly](nextflow.md).

## Prepare a sample sheet

TODO

## Initialize

Copy the configuration files to your project directory

```sh
champagne init --output /data/$USER/champagne_project
```

or if you do not use `--output`, your current working directory will be used as default:

```sh
champagne init
```

## Run

TODO preview, stub, mode=slurm

TODO required params

Run preview to view processes that will run:

```sh
champagne run --output /data/$USER/champagne_project -profile test -preview
```

Launch a stub run to view processes that will run and download containers:

```sh
champagne run --output /data/$USER/champagne_project -profile test,singularity -stub
```

Run the test dataset using the test profile:

```sh
champagne run --output /data/$USER/champagne_project -profile test,singularity
```

or explicitly specify the nextflow output directory and input:

```sh
champagne run --output /data/$USER/champagne_project -profile singularity --outdir results/test --input assets/samplesheet_test.csv
```

### Custom reference genome

TODO different required params

Create and use a custom reference genome:

```sh
champagne run --output /data/$USER/champagne_project -profile test -entry MAKE_REFERENCE
champagne run --output /data/$USER/champagne_project -profile test -c results/test/genome/custom_genome.config
```
