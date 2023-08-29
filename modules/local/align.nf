
process ALIGN_BLACKLIST {
    tag { sample_id }
    label 'align'

    input:
        tuple val(sample_id), path(fastq)
        path(blacklist_files)

    output:
        tuple val(sample_id), path("${sample_id}.no_blacklist.fastq.gz"), emit: reads

    script: // TODO use samtools -f4 for single-end and -f12 for paired to get unmapped reads https://broadinstitute.github.io/picard/explain-flags.html
    """
    bwa mem -t $task.cpus $params.align.blacklist $fastq |\
    samtools view -@ $task.cpus -f4 -b |\
    samtools bam2fq |\
    pigz -p $task.cpus > ${sample_id}.no_blacklist.fastq.gz
    """

    stub:
    """
    touch ${sample_id}.no_blacklist.fastq.gz
    """
}

process ALIGN_GENOME {
    tag { sample_id }
    label 'align'

    input:
        tuple val(sample_id), path(fastq)
        path(reference_files)

    output:
        tuple val(sample_id), path("${sample_id}.aligned.filtered.bam"), emit: bam

    script:
    """
    bwa mem -t ${task.cpus} ${params.align.genome} $fastq |\
    samtools sort -@ ${task.cpus} |\
    samtools view -b -q ${params.align.min_quality} -o ${sample_id}.aligned.filtered.bam
    """

    stub:
    """
    touch ${sample_id}.aligned.filtered.bam
    """

}

process INDEX_BAM {
    tag { sample_id }
    label 'align'

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("*.bam"), path("*.bai"), emit: bam
        path("*.idxstats"), emit: idxstats

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam $bam
    samtools index ${sample_id}.sorted.bam   # creates ${sample_id}.sorted.bam.bai
    samtools idxstats ${sample_id}.sorted.bam > ${sample_id}.sorted.bam.idxstats
    """

    stub:
    """
    touch ${sample_id}.sorted.bam ${sample_id}.sorted.bam.bai ${sample_id}.sorted.bam.idxstats
    """
}
