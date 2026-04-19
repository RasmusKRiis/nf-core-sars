
process IRMA {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'

    //conda "bioconda::irma=1.0.3"
    //container 'docker.io/rasmuskriis/cdc_irma_custom:1.0'
    container 'docker.io/cdcgov/irma:latest'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        //'https://depot.galaxyproject.org/singularity/irma:1.0.3--pl5321hdfd78af_0':
        //'biocontainers/irma:1.0.3--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)
   

    output:
    tuple val(meta), path("$meta.id/*.fasta") , emit: fasta
    tuple val(meta), path("$meta.id/*.bam"),  path("$meta.id/*.bai"), emit: bam
    tuple val(meta), path("$meta.id/*.vcf") , emit: vcf
    tuple val(meta), path("$meta.id/tables/READ_COUNTS.txt") , emit: read_count
    tuple val(meta), path("$meta.id/figures/*.pdf") , emit: figures
    tuple val(meta), path("$meta.id/amended_consensus/${meta.id}.fa") , emit: amended_consensus
    path("$meta.id/amended_consensus/${meta.id}.fa") , emit: amended_consensus_report
    tuple val(meta), path("$meta.id/secondary") , emit: secondary
  

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    IRMA CoV-minion-long-reads --external-config /project-bin/SARS-CoV-2-WGS-Nanopore.sh $fastq ${meta.id}
    """
}
