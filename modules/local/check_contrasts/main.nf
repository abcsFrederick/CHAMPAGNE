process CHECK_CONTRASTS {
    tag { contrasts }
    label 'process_single'

    container 'nciccbr/consensus_peaks:v1.1'

    input:
        path(samplesheet)
        path(contrasts)

    output:
        path("*.csv"),        emit: csv
        path("versions.yml"), emit: versions

    script:
    output_basename = "sample_contrasts"
    template 'check_contrasts.R'
}
