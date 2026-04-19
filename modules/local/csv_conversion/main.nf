
process CSV_CONVERSION {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'

    //conda "bioconda::irma=1.0.3"
    //container 'docker.io/rasmuskriis/cdc_irma_custom:1.0'
    container 'docker.io/rasmuskriis/nextclade-python'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        //'https://depot.galaxyproject.org/singularity/irma:1.0.3--pl5321hdfd78af_0':
        //'biocontainers/irma:1.0.3--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(csv)
   

    output:
    tuple val(meta), path("*nextclade_stats.csv") , emit: nextclade_stats
    tuple val(meta), path("*nextclade_mutations.csv") , emit: nextclade_mutations
    path("*nextclade_stats.csv") , emit: nextclade_stats_report
    path("*nextclade_mutations.csv") , emit: nextclade_mutations_report
    path("versions.yml"), emit: versions

  

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail


    #NEXTCLADE CONVERSION
    python3 /project-bin/csv_conversion_nextclade.py  \
                ${csv} \
                ${meta.id} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cat <<EOF > ${meta.id}_nextclade_stats.csv
    Sample,coverage,NC_Genome_MixedSites,NC_Genome_QC,NC_Genome_frameShifts,cdsCoverage
    ${meta.id},99,0,good,No frameShifts,S:1;ORF1a:1;ORF1b:1
    EOF

    cat <<EOF > ${meta.id}_nextclade_mutations.csv
    Sample,coverage,cdsCoverage,S_aaSubstitutions,ORF1a_aaSubstitutions,ORF1b_aaSubstitutions
    ${meta.id},99,S:1;ORF1a:1;ORF1b:1,N501Y,P3395H,P314L
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
