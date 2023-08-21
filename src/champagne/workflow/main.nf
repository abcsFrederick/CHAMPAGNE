log.info """\
         C H A M P A G N E 🍾
         ====================
         NF version   : $nextflow.version
         runName      : $workflow.runName
         username     : $workflow.userName
         configs      : $workflow.configFiles
         profile      : $workflow.profile
         cmd line     : $workflow.commandLine
         start time   : $workflow.start
         projectDir   : $workflow.projectDir
         launchDir    : $workflow.launchDir
         workdDir     : $workflow.workDir
         homeDir      : $workflow.homeDir
         reads        : ${params.reads}
         """
         .stripIndent()

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
