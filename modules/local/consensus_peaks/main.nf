process CONSENSUS_PEAKS {
    tag { meta.sample_basename }
    label 'peaks'
    label 'process_single'

    container "${params.containers.base}"

    input:
        tuple val(sample_tool), val(metas), val(peaks)
    output:
        tuple val(sample_tool), tuple val(metas), path("*.consensus_peaks.bed")

    script:
    // TODO assert that meta sample_basename is the same for all meta in metas
    """
    get_consensus_peaks.py --peakfiles ${peaks.join(' ')} --outbed ${sample_tool}.consensus_peaks.bed
    """

    stub:
    """
    touch ${sample_tool}.consensus_peaks.bed
    """
}
