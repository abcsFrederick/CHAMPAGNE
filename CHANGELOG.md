## CHAMPAGNE development version

### New features

- Implemented peak calling with sicer2, macs2, and gem. (#52)
- Added parameter options to skip QC, input normalization, and/or peak calling steps. (#72)

### Bug fixes

- CLI:
  - Error when biowulf-specific environment variables are not defined. (#54)
  - The host is now correctly detected as biowulf via `scontrol`. (#75)
- Containers:
  - When mount binding paths to containers, we use `--mount type=bind` for Singularity and `--volume` for Docker for compatibility. (#69)
  - Containers are now specified in process definitions instead of `withName`/`withLabel` for better control. (#69)
    - Shared containers are specified as parameters in the config file.

## CHAMPAGNE 0.1.0

### Quality control steps implemented for single-end reads

- Trim raw reads, FastQC on raw and trimmed reads, and FastQ Screen on trimmed reads.
- Exclude reads that align to blacklist regions, align remaining reads to the reference genome, and deduplicate.
- Preseq on aligned reads.
- Phantompeakqualtools on aligned and deduplicated reads.
- Process reads with deepTools: bam coverage to generate bigwigs for each sample, summarize all bigwigs, and compute matrices relative to TSSs and scaled to metagene regions.
- Generate plots with deepTools: PCA, profile, heatmap, spearman correlation, and fingerprint plots.
- Summarize all quality control steps in a MultiQC report.
- Input-normalize ChIP fragments for the next stage of the pipeline.
