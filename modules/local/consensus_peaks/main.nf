process CONSENSUS_PEAKS {
    tag "${meta.id}.${meta.group}"
    label 'peaks'
    label 'process_single'

    container "nciccbr/ccbr_ucsc:v1"

    input:
        tuple val(meta), val(peaks)
    output:
        tuple val(meta), path("*.consensus_peaks.bed"), emit: peaks

    script:
    if (peaks.size() > 1) {
        """
        get_consensus_peaks.py \\
            --peakfiles ${peaks.join(' ')} --outbed ${meta.id}.${meta.group}.consensus_peaks.bed
        """
    }
    // just copy the input if there's only one peak file
    else {
        """
        cp ${peaks.join(' ')} ${meta.id}.${meta.group}.consensus_peaks.bed
        """
    }

    stub:
    """
    touch ${meta.id}.${meta.group}.consensus_peaks.bed
    """
}
