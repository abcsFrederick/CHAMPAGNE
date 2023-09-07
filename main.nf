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
include { HANDLE_PRESEQ_ERROR      } from "./modules/local/qc.nf"
include { PARSE_PRESEQ_LOG         } from "./modules/local/qc.nf"
include { PHANTOM_PEAKS            } from "./modules/local/qc.nf"
include { PPQT_PROCESS             } from "./modules/local/qc.nf"
include { NGSQC_GEN                } from "./modules/local/qc.nf"
include { QC_STATS                 } from "./modules/local/qc.nf"
include { QC_TABLE                 } from "./modules/local/qc.nf"
include { MULTIQC                  } from "./modules/local/qc.nf"

include { ALIGN_BLACKLIST          } from "./modules/local/align.nf"
include { ALIGN_GENOME             } from "./modules/local/align.nf"
include { INDEX_BAM                } from "./modules/local/align.nf"

// TODO reorganize deeptools into separate subworkflow
include { BAM_COVERAGE             } from "./modules/local/deeptools.nf"
include { BIGWIG_SUM               } from "./modules/local/deeptools.nf"
include { BED_PROTEIN_CODING       } from "./modules/local/deeptools.nf"
include { COMPUTE_MATRIX           } from "./modules/local/deeptools.nf"
include { PLOT_FINGERPRINT         } from "./modules/local/deeptools.nf"
include { PLOT_CORRELATION         } from "./modules/local/deeptools.nf"
include { PLOT_PCA                 } from "./modules/local/deeptools.nf"
include { PLOT_HEATMAP             } from "./modules/local/deeptools.nf"
include { PLOT_PROFILE             } from "./modules/local/deeptools.nf"
include { NORMALIZE_INPUT          } from "./modules/local/deeptools.nf"

// TODO reorganize peak calling into separate subworkflow
include { CALC_GENOME_FRAC } from "./modules/local/peaks.nf"
include { SICER            } from "./modules/local/peaks.nf"

// MAIN WORKFLOW
workflow {
  INPUT_CHECK (
      file(params.input),
      params.seq_center
  )
  raw_fastqs = INPUT_CHECK.out.reads
  raw_fastqs.combine(Channel.value("raw")) | FASTQC_RAW
  raw_fastqs | TRIM_SE
  trimmed_fastqs = TRIM_SE.out
  trimmed_fastqs.combine(Channel.value("trimmed")) | FASTQC_TRIMMED
  trimmed_fastqs.combine(Channel.fromPath(params.fastq_screen.conf)) | FASTQ_SCREEN
  blacklist_files = Channel
                    .fromPath(params.align.blacklist_files)
                    .collect()
  ALIGN_BLACKLIST(trimmed_fastqs, blacklist_files)
  reference_files = Channel
                    .fromPath(params.align.reference_files)
                    .collect()
  ALIGN_GENOME(ALIGN_BLACKLIST.out.reads, reference_files)

  PRESEQ(ALIGN_GENOME.out.bam)
  // when preseq fails, write NAs for the stats that are calculated from its log
  PRESEQ.out.log
    .join(ALIGN_GENOME.out.bam, remainder: true)
    .branch { meta, preseq_log, bam_tuple ->
      failed: preseq_log == null
        return (tuple(meta, "nopresqlog"))
      succeeded: true
        return (tuple(meta, preseq_log))
    }.set{ preseq_logs }
  preseq_logs.failed | HANDLE_PRESEQ_ERROR
  preseq_logs.succeeded | PARSE_PRESEQ_LOG
  PARSE_PRESEQ_LOG.out.nrf
    .concat(HANDLE_PRESEQ_ERROR.out.nrf)
    .set{ preseq_nrf }

  chrom_sizes = Channel.fromPath(params.align.chrom_sizes)
  ALIGN_GENOME.out.bam.combine(chrom_sizes) | DEDUPLICATE
  DEDUPLICATE.out.bam | INDEX_BAM

  // NGSQC is seg faulting, see https://github.com/CCBR/CHAMPAGNE/issues/13
  //DEDUPLICATE.out.tag_align.combine(chrom_sizes) | NGSQC_GEN

  INDEX_BAM.out.bam | PHANTOM_PEAKS
  PPQT_PROCESS(PHANTOM_PEAKS.out.fraglen)
  frag_lengths = PPQT_PROCESS.out.fraglen
  QC_STATS(
    raw_fastqs,
    ALIGN_GENOME.out.flagstat,
    DEDUPLICATE.out.flagstat,
    preseq_nrf,
    PHANTOM_PEAKS.out.spp,
    frag_lengths
  )
  QC_TABLE(QC_STATS.out.collect())

  // Deeptools
  BAM_COVERAGE(INDEX_BAM.out.bam, frag_lengths)
  BIGWIG_SUM(BAM_COVERAGE.out.bigwig.collect())
  BIGWIG_SUM.out.array.combine(Channel.from('heatmap', 'scatterplot')) | PLOT_CORRELATION
  BIGWIG_SUM.out.array | PLOT_PCA

  // Create channel: [ meta, [ ip_bam, control_bam ] [ ip_bai, control_bai ] ]
  ch_genome_bam_bai = INDEX_BAM.out.bam
  ch_genome_bam_bai
      .combine(ch_genome_bam_bai)
      .map {
          meta1, bam1, bai1, meta2, bam2, bai2 ->
              meta1.control == meta2.id ? [ meta1, [ bam1, bam2 ], [ bai1, bai2 ] ] : null
      }
      .set { ch_ip_ctrl_bam_bai }

  PLOT_FINGERPRINT(ch_ip_ctrl_bam_bai)
  BED_PROTEIN_CODING(Channel.fromPath(params.gene_info))
  COMPUTE_MATRIX(BAM_COVERAGE.out.bigwig.collect(),
                 BED_PROTEIN_CODING.out.bed.combine(Channel.from('metagene','TSS'))
  )
  PLOT_HEATMAP(COMPUTE_MATRIX.out.mat)
  PLOT_PROFILE(COMPUTE_MATRIX.out.mat)

  // Create channel: [ meta, ip_bw, control_bw ]
  BAM_COVERAGE.out.meta
      .merge(BAM_COVERAGE.out.bigwig)
      .set { bigwigs }
  bigwigs
      .combine(bigwigs)
      .map {
        meta1, bw1, meta2, bw2 ->
            meta1.control == meta2.id ? [ meta1, bw1, bw2 ] : null
      }
      .set { ch_ip_ctrl_bigwig }

  MULTIQC(
    Channel.fromPath(params.multiqc_config),
    FASTQC_RAW.out.zip.collect(),
    FASTQC_TRIMMED.out.zip.collect(),
    FASTQ_SCREEN.out.screen.collect(),
    //NGSQC_GEN
    DEDUPLICATE.out.flagstat.collect(),
    PHANTOM_PEAKS.out.spp.collect(),
    QC_TABLE.out.txt,
    PLOT_FINGERPRINT.out.matrix.collect(),
    PLOT_FINGERPRINT.out.metrics.collect(),
    PLOT_CORRELATION.out.tab.collect(),
    PLOT_PCA.out.tab.collect(),
    PLOT_PROFILE.out.tab.collect()
  )

  NORMALIZE_INPUT(ch_ip_ctrl_bigwig)

  // peak calling

  genome_frac = CALC_GENOME_FRAC(chrom_sizes)
  // create channel with [ meta, chip_tag, input_tag, fraglen, genome_frac]
  DEDUPLICATE.out.tag_align
    .combine(DEDUPLICATE.out.tag_align)
    .map {
        meta1, tag1, meta2, tag2 ->
            meta1.control == meta2.id ? [ meta1, tag1, tag2 ]: null
    }
    .join(frag_lengths)
    .combine(genome_frac)
    .set { ch_ip_ctrl_tagalign }

  ch_ip_ctrl_tagalign | SICER

}
