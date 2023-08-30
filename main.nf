log.info """\
         CHAMPAGNE 🍾
         =============
         NF version   : $nextflow.version
         runName      : $workflow.runName
         username     : $workflow.userName
         configs      : $workflow.configFiles
         profile      : $workflow.profile
         cmd line     : $workflow.commandLine
         start time   : $workflow.start
         projectDir   : $workflow.projectDir
         launchDir    : $workflow.launchDir
         workDir      : $workflow.workDir
         homeDir      : $workflow.homeDir
         input        : ${params.input}
         """
         .stripIndent()

// SUBWORKFLOWS
include { INPUT_CHECK              } from './subworkflows/local/input_check.nf'

// MODULES
include { TRIM_SE                  } from "./modules/local/trim.nf"
// TODO reorganize QC into separate subworkflow
include { FASTQC as FASTQC_RAW     } from "./modules/local/qc.nf"
include { FASTQC as FASTQC_TRIMMED } from "./modules/local/qc.nf"
include { FASTQ_SCREEN             } from "./modules/local/qc.nf"
include { DEDUPLICATE              } from "./modules/local/qc.nf"
include { PRESEQ                   } from "./modules/local/qc.nf"
include { PHANTOM_PEAKS            } from "./modules/local/qc.nf"
include { NGSQC_GEN                } from "./modules/local/qc.nf"

include { ALIGN_BLACKLIST          } from "./modules/local/align.nf"
include { ALIGN_GENOME             } from "./modules/local/align.nf"
include { INDEX_BAM                } from "./modules/local/align.nf"

// TODO reorganize deeptools into separate subworkflow
include { BAM_COVERAGE             } from "./modules/local/deeptools.nf"
include { BIGWIG_SUM               } from "./modules/local/deeptools.nf"
include { ARRAY_PLOTS              } from "./modules/local/deeptools.nf"
include { FINGERPRINT              } from "./modules/local/deeptools.nf"
include { BED_PROTEIN_CODING       } from "./modules/local/deeptools.nf"
include { COMPUTE_MATRIX           } from "./modules/local/deeptools.nf"
include { PLOT_HEATMAP             } from "./modules/local/deeptools.nf"
include { PLOT_PROFILE             } from "./modules/local/deeptools.nf"

// MAIN WORKFLOW
workflow {

  INPUT_CHECK (
      file(params.input),
      params.seq_center
  )
  raw_fastqs = INPUT_CHECK.out.reads
  raw_fastqs | TRIM_SE
  trimmed_fastqs = TRIM_SE.out
  trimmed_fastqs.combine(Channel.value("trimmed")) | FASTQC_TRIMMED
  //trimmed_fastqs = TRIM_SE.out | map { it -> it[1] } | collect // obtain list of just fastqs without sample_ids
  trimmed_fastqs.combine(Channel.fromPath(params.fastq_screen.conf)) | FASTQ_SCREEN
  blacklist_files = Channel
                    .fromPath("${params.align.index_dir}${params.align.blacklist}*")
                    .collect()
  ALIGN_BLACKLIST(trimmed_fastqs, blacklist_files)
  reference_files = Channel
                    .fromPath("${params.align.index_dir}${params.align.genome}*")
                    .collect()
  ALIGN_GENOME(ALIGN_BLACKLIST.out, reference_files)
  //PRESEQ(ALIGN_GENOME.out.bam)
  chrom_sizes = Channel.fromPath("${params.align.index_dir}${params.align.chrom_sizes}")
  ALIGN_GENOME.out.bam.combine(chrom_sizes) | DEDUPLICATE
  DEDUPLICATE.out.bam | INDEX_BAM
  //DEDUPLICATE.out.tag_align.combine(chrom_sizes) | NGSQC_GEN
  INDEX_BAM.out.bam | PHANTOM_PEAKS
  BAM_COVERAGE(INDEX_BAM.out.bam, PHANTOM_PEAKS.out.ppqt)
  BIGWIG_SUM(BAM_COVERAGE.out.bigwig.collect())
  BIGWIG_SUM.out.array | ARRAY_PLOTS

  // Create channels: [ meta, [ ip_bam, control_bam ] [ ip_bai, control_bai ] ]
  ch_genome_bam_bai = INDEX_BAM.out.bam
  ch_genome_bam_bai
      .combine(ch_genome_bam_bai)
      .map {
          meta1, bam1, bai1, meta2, bam2, bai2 ->
              meta1.control == meta2.id ? [ meta1, [ bam1, bam2 ], [ bai1, bai2 ] ] : null
      }
      .set { ch_ip_control_bam_bai }

  //FINGERPRINT(ch_ip_control_bam_bai) // TODO https://github.com/CCBR/Dockers/issues/12
  BED_PROTEIN_CODING(Channel.fromPath(params.gene_info))
  COMPUTE_MATRIX(BAM_COVERAGE.out.bigwig.collect(),
                 BED_PROTEIN_CODING.out.bed.combine(Channel.from('metagene','TSS'))
  )
  PLOT_HEATMAP(COMPUTE_MATRIX.out)
  PLOT_PROFILE(COMPUTE_MATRIX.out)
}
