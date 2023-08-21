
include { TRIM_SE } from "./modules/local/trim.nf"

workflow convert2fasta {
  Channel.fromPath(params.input) | any2fasta | view
}

workflow {
  raw_fastqs = Channel
                    .fromPath(params.reads)
                    .map { file -> tuple(file.simpleName, file) }
  TRIM_SE(raw_fastqs)
}
