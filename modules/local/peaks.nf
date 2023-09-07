
process CALC_GENOME_FRAC {
    label 'peaks'

    input:
        path(chrom_sizes)

    output:
        env(genome_frac), emit: genome_frac

    script:
    """
    genome_frac=`calc_effective_genome_fraction.py ${params.align.effective_genome_size} ${chrom_sizes} ${params.deeptools.excluded_chroms}`
    echo \$genome_frac
    """

    stub:
    """
    genome_frac=0.75
    echo \$genome_frac
    """
}

process SICER {
    label 'peaks'

    input:
        tuple val(meta), path(chip), path(input), val(fraglen), val(genome_frac)

    //output:

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
    # TODO
    """
}
/*
process SICER_NO_CTRL { // TODO just use one process with optional ctrl input
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
