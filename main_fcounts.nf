nextflow.enable.dsl = 2

// Crear timestamp para versionar resultados
def timestamp = new Date().format("yyyyMMdd_HHmmss")

process fastqc {
    publishDir "results/fastqc_${timestamp}", mode:'copy', overwrite: false

    input:
        tuple val(sample), path(forward), path(reverse)

    output:
        path "*_fastqc.*"

    script:
    """
    fastqc ${forward} ${reverse} --outdir ./
    """
}

process fastqc_trimmed {
    publishDir "results/fastqc_posttrim_${timestamp}", mode:'copy', overwrite: false

    input:
        tuple val(sample), path(forward_trimmed), path(reverse_trimmed), path(unused1), path(unused2)

    output:
        path "*_fastqc.*"

    script:
    """
    fastqc ${forward_trimmed} ${reverse_trimmed} --outdir ./
    """
}

process trimmomatic {
    publishDir "results/trimmomatic_${timestamp}", mode:'copy', overwrite: false

    input:
        tuple val(sample), path(forward), path(reverse)
        path adapters_file

    output:
        tuple val(sample), \
              path("${sample}_paired_forward.fastq"), path("${sample}_unpaired_forward.fastq"), \
              path("${sample}_paired_reverse.fastq"), path("${sample}_unpaired_reverse.fastq")

    script:
    """
    trimmomatic PE -threads ${task.cpus} \
        ${forward} ${reverse} \
        ${sample}_paired_forward.fastq ${sample}_unpaired_forward.fastq \
        ${sample}_paired_reverse.fastq ${sample}_unpaired_reverse.fastq \
        ILLUMINACLIP:${adapters_file}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process hisat2 {
    label 'high_memory'
    publishDir "results/hisat2_${timestamp}", mode: 'copy', overwrite: false

    input:
        tuple val(sample), path(forward_paired), path(unused_forward), path(reverse_paired), path(unused_reverse)
        path index_files

    output:
        tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai")

    script:
    """
    hisat2 -p ${task.cpus} \
           -x genome \
           -1 ${forward_paired} \
           -2 ${reverse_paired} | \
        samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam

    samtools index ${sample}.sorted.bam
    """
}

process featurecounts {
    conda 'envs/featurecounts.yml'
    publishDir "results/featurecounts_${timestamp}", mode: 'copy', overwrite: false

    input:
        tuple val(sample), path(bam), path(bai)
        path gtf_annotation

    output:
        tuple val(sample), path("${sample}_counts.txt")

    script:
    """
    featureCounts -T ${task.cpus} \
        -a ${gtf_annotation} \
        -o ${sample}_counts.txt \
        -p -B -C \
        ${bam}
    """
}

workflow {
    adapters = file("TruSeq3-PE.fa")
    hisat2_index = file("reference/grch38/*.ht2")
    gtf_annotation = file("reference/Homo_sapiens.GRCh38.109.gtf")

    Channel.fromFilePairs("data/fastq/*_{forward,reverse}.fastq", flat: true)
           .set { paired_reads }

    fastqc(paired_reads)

    trimmed_reads = trimmomatic(paired_reads, adapters)

    fastqc_trimmed(trimmed_reads)

    aligned = hisat2(trimmed_reads, hisat2_index)

    quantified = featurecounts(aligned, gtf_annotation)
}

