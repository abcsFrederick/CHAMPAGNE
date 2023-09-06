
process calc_effective_genome_fraction {

}

process SICER {
    label 'peaks'

    input:
        tuple val(meta), path(chip), path(input), val(fraglen), val(genome_frac)

    script:
    """
    sicer \\
      -t ${chip} \\
      -c ${input} \\
      -s ${params.genome} \\
      -rt 100 \\
      -w 300 \\
      -f ${fraglen} \\
      -egf ${genome_frac} \\ // TODO
      -g 600 \\
      -fdr 1E-2 \\
      -cpu ${task.cpus}
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
