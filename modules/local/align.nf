
process ALIGN_BLACKLIST {
    tag { meta.id }
    label 'align'

    input:
        tuple val(meta), path(fastq)
        path(blacklist_files)

    output:
        tuple val(meta), path("${meta.id}.no_blacklist.fastq.gz"), emit: reads

    script: // TODO use samtools -f4 for single-end and -f12 for paired to get unmapped reads https://broadinstitute.github.io/picard/explain-flags.html
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bwa mem -t $task.cpus $params.align.blacklist $fastq |\
    samtools view -@ $task.cpus -f4 -b |\
    samtools bam2fq |\
    pigz -p $task.cpus > ${prefix}.no_blacklist.fastq.gz
    """

    stub:
    """
    touch ${meta.id}.no_blacklist.fastq.gz
    """
}

process ALIGN_GENOME {
    tag { meta.id }
    label 'align'

    input:
        tuple val(meta), path(fastq)
        path(reference_files)

    output:
        tuple val(meta), path("*.aligned.filtered.bam"), emit: bam
        path("*.aligned.filtered.bam.flagstat"), emit: flagstat

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bwa mem -t ${task.cpus} ${params.align.genome} $fastq |\
    samtools sort -@ ${task.cpus} |\
    samtools view -b -q ${params.align.min_quality} -o ${prefix}.aligned.filtered.bam
    samtools flagstat ${prefix}.aligned.filtered.bam > ${prefix}.aligned.filtered.bam.flagstat
    """

    stub:
    """
    touch ${meta.id}.aligned.filtered.bam
    """

}

process INDEX_BAM {
    tag { meta.id }
    label 'align'

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
