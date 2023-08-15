# 4. Expected Outputs

TODO The following directories are created under the WORKDIR/results directory:

- alignment_stats: this directory include information on the alignment of each sample
- peaks: this directory contains a sub-directory that relates to the quality threshold used.
  - quality threshold
    - contrasts: this directory includes the contrasts for each line listed in the contrast manifest
    - peak_caller: this directory includes all peak calls from each peak_caller (SEACR, MACS2, GOPEAKS) for each sample
      - annotation
        - go_enrichment: this directory includes gene set enrichment pathway predictions
        - homer: this directory includes the annotation output from HOMER
        - rose: this directory includes the annotation output from ROSE

```
# TODO use tree command to show working directory structure
```
