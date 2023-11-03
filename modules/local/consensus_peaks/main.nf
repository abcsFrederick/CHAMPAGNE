process CONSENSUS_PEAKS {
    tag "${sample}_${tool}"
    label 'peaks'
    label 'process_single'

    container "nciccbr/ccbr_ucsc:v1"

    input:
        tuple val(sample), val(tool), val(metas), val(peaks)
    output:
        tuple val(sample), val(tool), path("*.consensus_peaks.bed"), emit: peaks

    script:
    // assert that meta sample_basename is the same for all
    assert metas.collect{ it.sample_basename }.toSet().size() == 1
    if (metas.size() > 1) {
        """
        get_consensus_peaks.py --peakfiles ${peaks.join(' ')} --outbed ${sample}.${tool}.consensus_peaks.bed
        """
    }
    // just copy the input if there's only one replicate
    else {
        """
        cp ${peaks} ${sample}.${tool}.consensus_peaks.bed
        """
    }

    stub:
    """
    touch ${sample}.${tool}.consensus_peaks.bed
    """
}
