
process CALC_GENOME_FRAC {
    label 'peaks'

    container = "${params.containers.base}"

    input:
        path(chrom_sizes)

    output:
        env(genome_frac), emit: genome_frac

    script:
    """
    genome_frac=`calc_effective_genome_fraction.py ${params.genomes[ params.genome ].effective_genome_size} ${chrom_sizes} ${params.deeptools.excluded_chroms}`
    echo \$genome_frac
    """

    stub:
    """
    genome_frac=0.75
    echo \$genome_frac
    """
}


process MACS_BROAD {
    tag { meta.id }
    label 'peaks'

    container = "${params.containers.macs2}"

    input:
        tuple val(meta), path(chip), path(input), val(fraglen), val(genome_frac)

    output:
        tuple val(meta), path("${meta.id}_peaks.broadPeak"), val("${task.process.tokenize(':')[-1].toLowerCase()}"), emit: peak
        path("${meta.id}_peaks.xls")
        path("${meta.id}_peaks.gappedPeak")

    script:
    """
    macs2 callpeak \\
      -t ${chip} \\
      -c ${input} \\
      -g ${params.genomes[ params.genome ].effective_genome_size} \\
      -n ${meta.id} \\
      --extsize ${fraglen} \\
      --nomodel \\
      -q 0.01 \\
      --keep-dup='all' \\
      --broad \\
      --broad-cutoff 0.01
    """

    stub:
    """
    for ext in xls broadPeak gappedPeak; do
        touch ${meta.id}_peaks.\${ext}
    done
    """

}

process MACS_NARROW {
    tag { meta.id }
    label 'peaks'

    container = "${params.containers.macs2}"

    input:
        tuple val(meta), path(chip), path(input), val(fraglen), val(genome_frac)

    output:
        tuple val(meta), path("${meta.id}_peaks.narrowPeak"), val("${task.process.tokenize(':')[-1].toLowerCase()}"), emit: peak
        path("${meta.id}_peaks.xls")
        path("${meta.id}_summits.bed")

    script:
    """
    macs2 callpeak \\
      -t ${chip} \\
      -c ${input} \\
      -g ${params.genomes[ params.genome ].effective_genome_size} \\
      -n ${meta.id} \\
      --extsize ${fraglen} \\
      --nomodel \\
      -q 0.01 \\
      --keep-dup='all'
    """

    stub:
    """
    for ext in peaks.xls peaks.narrowPeak summits.bed; do
        touch ${meta.id}_\${ext}
    done
    """
}

process SICER {
    tag { meta.id }
    label 'peaks'
    label 'process_high'

    container = "${params.containers.sicer}"

    input:
        tuple val(meta), path(chip), path(input), val(fraglen), val(genome_frac)

    output:
        tuple val(meta), path("*island.bed"), val("${task.process.tokenize(':')[-1].toLowerCase()}"), emit: peak
        //tuple path("*.scoreisland"), path("*normalized.wig"), path("*islands-summary")

    script:
    """
    sicer \\
      -t ${chip} \\
      -c ${input} \\
      -s ${params.genome} \\
      -rt 100 \\
      -w 300 \\
      -f ${fraglen} \\
      -egf ${genome_frac} \\
      -g 600 \\
      -fdr 1E-2 \\
      -cpu ${task.cpus}
    """

    stub:
    """
    for ext in scoreisland normalized.wig islands-summary island.bed; do
        touch ${chip.baseName}.\${ext}
    done
    """
}

process CONVERT_SICER { // https://github.com/CCBR/Pipeliner/blob/86c6ccaa3d58381a0ffd696bbf9c047e4f991f9e/Rules/ChIPseq.snakefile#L389-L431
    tag { meta.id }
    label 'peaks'
    label 'process_single'

    container = "${params.containers.base}"

    input:
        tuple val(meta), path(sicer_peaks), val(peak_tool)

    output:
        tuple val(meta), path("${sicer_peaks.baseName}.converted.bed"), val(peak_tool), emit: peak

    script:
    $/
    #!/usr/bin/env python
    import math
    with open("${sicer_peaks}",'r') as f:
        intxt = f.readlines()
    # input columns if input-normalized: chrom, start, end, ChIP tag count, control tag count, p-value, fold-enrichment, q-value
    # input columns if no input: chrom, start, end, score
    outBroadPeak = [None] * len(intxt)
    outBed = [None] * len(intxt)
    for i in range(len(intxt)):
        tmp = intxt[i].strip().split('\t')
        if len(tmp) == 8:
            # assuming that a p-value/q-value of 0 is super significant, -log10(1e-500)
            if tmp[5] == "0.0":
                pval="500"
            else:
                pval = str(-(math.log10(float(tmp[5]))))
            if tmp[7] == "0.0":
                qval="500"
                qvalScore="5000"
            else:
                qval = str(-(math.log10(float(tmp[7]))))
                qvalScore = str(int(-10*math.log10(float(tmp[7]))))
            outBroadPeak[i] = "\t".join(tmp[0:3] + ["Peak"+str(i+1),tmp[3],".", tmp[6], pval, qval])
            outBed[i] = "\t".join(tmp[0:3] + ["Peak"+str(i+1),qvalScore])
        else:
            score = str(int(float(tmp[3])))
            outBed[i] = "\t".join(tmp[0:3] + ["Peak"+str(i+1),score])
    # output broadPeak columns: chrom, start, end, name, ChIP tag count, strand, fold-enrichment, -log10 p-value, -log10 q-value
    with open("${sicer_peaks.baseName}.converted.bed",'w') as g:
        g.write( "\n".join(outBed) )
    if outBroadPeak[0] != None:
        with open("${sicer_peaks.baseName}.converted_sicer.broadPeak", 'w') as h:
            h.write( "\n".join(outBroadPeak) )
    /$

    stub:
    """
    touch ${sicer_peaks.baseName}.converted.bed ${sicer_peaks.baseName}.converted_sicer.broadPeak
    """

}
/*
process SICER_NO_CTRL { // TODO just use one process for each peak caller with optional ctrl input
    script:
    """
    sicer \\
      -t ${chip} \\
      -s ${params.genome} \\
      -rt 100 \\
      -w 300 \\
      -f ${fraglen} \\
      -egf ${genome_frac} \\
      -g 600 \\
      -e 100 \\
      -cpu ${task.cpus}
    """
}
*/
process GEM {
    tag { meta.id }
    label 'peaks'
    label 'process_high'

    container = "${params.containers.gem}"

    input:
        tuple val(meta), path(chip), path(input), path(read_dists), path(chrom_sizes), path(chrom_dir)

    output:
        tuple val(meta), path("${meta.id}/*GEM_events.narrowPeak"), val("${task.process.tokenize(':')[-1].toLowerCase()}"), emit: peak
        tuple val(meta), path("${meta.id}/*GEM_events.txt"), val("${task.process.tokenize(':')[-1].toLowerCase()}"), emit: event

    script:
    // $GEMJAR is defined in the docker container
    """
    java -Xmx30g -jar \$GEMJAR \\
      --t ${task.cpus} \\
      --d ${read_dists} \\
      --g ${chrom_sizes} \\
      --genome ${chrom_dir} \\
      --expt ${chip} \\
      --ctrl ${input} \\
      --out ${meta.id} \\
      --k_min 6 \\
      --k_max 13 \\
      --outNP \\
      --nrf
    """

    stub:
    """
    mkdir ${meta.id}
    for ext in txt narrowPeak; do
        touch ${meta.id}/${meta.id}.GEM_events.\$ext
    done
    """
}

process FRACTION_IN_PEAKS {
    tag { meta.id }
    label 'peaks'
    label 'process_single'

    container "${params.containers.frip}"

    input:
        tuple val(meta), path(dedup_bam), path(dedup_bai), path(peaks), val(peak_tool), path(chrom_sizes)

    output:
        path("*.FRiP.txt")

    script:
    """
    export OPENBLAS_NUM_THREADS='1'
    frip.py -p ${peaks} -b ${dedup_bam} -g ${chrom_sizes} -t ${peak_tool} -s ${meta.id} -o ${meta.id}_${peak_tool}.FRiP.txt
    """

    stub:
    """
    touch ${meta.id}_${peak_tool}.FRiP.txt
    """
}

process CONCAT_FRIPS {
    label 'peaks'
    label 'process_single'

    container "${params.containers.base}"

    input:
        path(frips)

    output:
        path("FRiP.txt")

    script:
    """
    header=`head -n 1 ${frips.get(1)}`
    echo \$header > FRiP.txt
    for f in ${frips}; do
        tail -n 1 \$f >> FRiP.txt
    done
    """
    stub:
    """
    touch FRiP.txt
    """

}

process PLOT_FRIP {
    tag { meta.id }
    label 'peaks'
    label 'process_single'

    container "${params.containers.r}"

    input:
        path(frips)

    output:
        path("*.png")

    script:
    """
    plot_frip.R ${frips}
    """

    stub:
    """
    touch blank.png
    """
}
