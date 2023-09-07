## CHAMPAGNE development version

### New features

- New GitHub Actions workflow to launch a stub run from the champagne CLI.

### Bug fixes

- CLI error when biowulf-specific environment variables are not defined. (#54)

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
