
process FILTER_BLACKLIST { // TODO: refactor this as a subworkflow that combines alignment and filtering processes
    tag { meta.id }
    label 'align'
    label 'process_higher'

    container = "${params.containers.base}"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("${meta.id}.no_blacklist.bam"), emit: bam

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def filter_flag = meta.single_end ? '4' : '12'
    """
    samtools view \\
      -@ ${task.cpus} \\
      -f${filter_flag} \\
      -b \\
      ${prefix}.bam > ${meta.id}.no_blacklist.bam
    """

    stub:
    """
    touch ${meta.id}.no_blacklist.bam
    """
}

process BAM_TO_FASTQ {
    tag { meta.id }
    label 'align'
    container "${ meta.single_end ? params.containers.base : params.containers.picard }"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*.R?.fastq*"), emit: reads

    script:
    if (meta.single_end) {
        """
        samtools bam2fq ${bam} | pigz -p ${task.cpus} > ${bam.baseName}.R1.fastq.gz
        """
    } else {
        """
        picard -Xmx${task.memory.toGiga()}G SamToFastq \\
            --VALIDATION_STRINGENCY SILENT \\
            --INPUT ${bam} \\
            --FASTQ ${bam.baseName}.R1.fastq \\
            --SECOND_END_FASTQ ${bam.baseName}.R2.fastq \\
            --UNPAIRED_FASTQ ${bam.baseName}.unpaired.fastq
        pigz -p ${task.cpus} *.fastq
        """
    }
    stub:
    """
    touch ${bam.baseName}.R1.fastq
    """
}

process ALIGN_GENOME { // TODO refactor as subworkflow
    tag { meta.id }
    label 'align'
    label 'process_higher'

    container = "${params.containers.base}"

    input:
        tuple val(meta), path(fastq)
        tuple val(meta_ref), path(reference_files)

    output:
        tuple val(meta), path("*.aligned.filtered.bam"), emit: bam
        path("*.aligned.filtered.bam.flagstat"), emit: flagstat

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def filter = ''
    if (meta.single_end) {
        filter = """samtools view \\
            -@ ${task.cpus} \\
            -q ${params.align_min_quality} \\
            -b \\
            ${prefix}.sorted.bam > ${prefix}.aligned.filtered.bam
        """
    } else {
        filter = """bam_filter_by_mapq.py \\
            -q ${params.align_min_quality} \\
            -i ${prefix}.sorted.bam \\
            -o ${prefix}.aligned.filtered.bam
        """
    }
    """
    # current working directory is a tmpdir when 'scratch' is set
    TMP=tmp/
    mkdir \$TMP
    trap 'rm -rf "\$TMP"' EXIT

    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa mem -t ${task.cpus} \$INDEX ${fastq} > ${prefix}.bam
    samtools sort \\
      -@ ${task.cpus} \\
      -m 2G \\
      -T \$TMP \\
      ${prefix}.bam > ${prefix}.sorted.bam
    samtools index ${prefix}.sorted.bam
    ${filter}
    samtools flagstat ${prefix}.aligned.filtered.bam > ${prefix}.aligned.filtered.bam.flagstat
    """

    stub:
    """
    touch ${meta.id}.aligned.filtered.bam ${meta.id}.aligned.filtered.bam.flagstat
    """

}

process INDEX_BAM {
    tag { meta.id }
    label 'align'
    label 'process_higher'

    container = "${params.containers.base}"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*.bam"), path("*.bai"), emit: bam
        path("*.idxstats"), emit: idxstats

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam $bam
    samtools index ${prefix}.sorted.bam   # creates ${prefix}.sorted.bam.bai
    samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
    """

    stub:
    """
    for ext in sorted.bam sorted.bam.bai sorted.bam.idxstats; do
        touch ${meta.id}.\$ext
    done
    """
}
