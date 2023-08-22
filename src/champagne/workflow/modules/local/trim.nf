// https://github.com/CCBR/ChiP-seek/blob/9ba449e4855f9710e86f2db7c1d9560de634b3f1/workflow/rules/align.smk#L21
// https://github.com/nf-core/ampliseq/blob/dev/subworkflows/local/cutadapt_workflow.nf
process TRIM_SE {
  tag { sample_id }
  publishDir "$params.outdir/$sample_id", mode: 'copy'

  input:
    tuple val(sample_id), path(fastq)

  output:
    tuple val(sample_id), path("${sample_id}.cutadapt.fastq")

  script:
  """
  cutadapt \
    --nextseq-trim=2 \
    --trim-n \
    -n 5 -O 5 \
    -q ${params.leadingquality},${params.trailingquality} \
    -m ${params.minlen} \
    -b file:${params.adapters} \
    -j $task.cpus \
    -o ${sample_id}.cutadapt.fastq \
    $fastq
  """
}
