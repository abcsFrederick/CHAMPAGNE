# nextflow pipeline

You can run the nextflow pipeline directly by specifying this GitHub repo.
You will need nextflow and either singularity or docker installed.
In this case you don't need to run `champagne init` first,
as the config files will be accessed directly from the GitHub repo.

```sh
nextflow run CCBR/CHAMPAGNE -profile test,singularity
```

You can specify a specific version, tag, or branch with `-r`:

```sh
nextflow run CCBR/CHAMPAGNE -r v1.0.0 -profile test,singularity
```

Create and use a custom reference genome:

```sh
nextflow run CCBR/CHAMPAGNE -profile test -entry MAKE_REFERENCE
nextflow run CCBR/CHAMPAGNE -profile test -c results/test/genome/custom_genome.config
```
