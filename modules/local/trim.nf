// https://github.com/CCBR/ChiP-seek/blob/9ba449e4855f9710e86f2db7c1d9560de634b3f1/workflow/rules/align.smk#L21
// https://github.com/nf-core/ampliseq/blob/dev/subworkflows/local/cutadapt_workflow.nf
process TRIM_SE {
  tag { sample_id }
  publishDir "$params.outdir/trimmed", mode: 'copy'

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
    -q ${params.cutadapt.leadingquality},${params.cutadapt.trailingquality} \
    -m ${params.cutadapt.minlen} \
    -b file:${params.cutadapt.adapters} \
    -j $task.cpus \
    -o ${sample_id}.cutadapt.fastq \
    $fastq
  """

  stub:
  """
  touch ${sample_id}.cutadapt.fastq
  """
}
