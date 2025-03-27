nextflow.enable.dsl = 2

process fastqc {
    publishDir "results/fastqc", mode:'copy'

    input:
        tuple val(sample), path(forward), path(reverse)

    output:
        path "*_fastqc.*"

    script:
    """
    fastqc ${forward} ${reverse} --outdir ./
    """
}

process trimmomatic {
    publishDir "results/trimmomatic", mode:'copy'

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
    publishDir "results/hisat2", mode: 'copy'

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

process stringtie {
    publishDir "results/stringtie", mode: 'copy'

    input:
        tuple val(sample), path(bam), path(bai)
        path gtf_annotation

    output:
        tuple val(sample), path("${sample}.gtf")

    script:
    """
    stringtie ${bam} \
        -G ${gtf_annotation} \
        -o ${sample}.gtf \
        -p ${task.cpus} \
        -e -B
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

    aligned = hisat2(trimmed_reads, hisat2_index)
    stringtie(aligned, gtf_annotation)
}

