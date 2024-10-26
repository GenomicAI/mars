#!/usr/bin/env nextflow
// Author : Shanika Amarasoma, Nuzla Ismail
// Date : October 27, 2024
// Run in singularity container
// singularity exec -e -B $PWD mars_latest.sif nextflow run example.nf

nextflow.enable.dsl=2

// Define input files and parameters
params.ref = "${baseDir}/NC_000020.11.fa"
params.sample = "${baseDir}/HG00096.fa"
params.prefix = "HG00096"
params.num_reads = 300000
params.length = 100
params.outdir = "./"

process SIMULATE_READS {
    input:
    path params.sample
    output:
    path "${params.outdir}/${params.prefix}_reads_R1.fq", emit: read1
    path "${params.outdir}/${params.prefix}_reads_R2.fq", emit: read2

    script:
    """
    wgsim -r 0 -e 0.001 -N ${params.num_reads} -1 ${params.length} -2 ${params.length} ${params.sample} \
        "${params.outdir}/${params.prefix}_reads_R1.fq" "${params.outdir}/${params.prefix}_reads_R2.fq"
    """
}

process ALIGN_READS {
    input:
    path read1
    path read2
    path params.ref
    output:
    path "${params.outdir}/${params.prefix}.bam"

    script:
    """
    bwa mem -R "@RG\\tID:${params.prefix}\\tSM:${params.prefix}\\tLB:L1" ${params.ref} ${read1} ${read2} | \
        samtools sort -o "${params.outdir}/${params.prefix}.bam"
    """
}

process CALL_VARIANTS {
    input:
    path "${params.outdir}/${params.prefix}.bam"
    path params.ref
    output:
    path "${params.outdir}/${params.prefix}.vcf.gz"

    script:
    """
    samtools index "${params.outdir}/${params.prefix}.bam"
    bcftools mpileup -f ${params.ref} "${params.outdir}/${params.prefix}.bam" | \
        bcftools call -mv -Oz -o "${params.outdir}/${params.prefix}.vcf.gz"
    bcftools index "${params.outdir}/${params.prefix}.vcf.gz"
    """
}

// Workflow definition
workflow {
    read_simulation = SIMULATE_READS(params.sample)
    aligned_reads = ALIGN_READS(read_simulation.read1, read_simulation.read2, params.ref)
    CALL_VARIANTS(aligned_reads, params.ref)
}
