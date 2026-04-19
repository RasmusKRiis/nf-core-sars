
process TABLELOOKUP {
    tag "$meta.id"
    label 'process_single'
    errorStrategy 'ignore'
    
    

    //conda "bioconda::blast=2.15.0"
    container 'docker.io/rasmuskriis/blast_python_pandas:amd64'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory

    input:
    tuple val(meta), path(csv)
    path(spike)
    path(rdrp)
    path(clpro)
   

    output:
    //tuple val(meta), path("*.txt"), emit: genotype
    //tuple val(meta), path("*.csv"), emit: genotype_file
    tuple val(meta), path("*.csv"), emit: resistance_mutations
    path("*.csv"), emit: resistance_mutations_report
    path("versions.yml"), emit: versions



    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    #Resistance mutation lookup
    python3 /project-bin/lookup_mutations.py  \
                ${csv} \
                ${meta.id} \
                ${spike} \
                ${rdrp} \
                ${clpro} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cat <<EOF > ${meta.id}_resistance_mutations.csv
    Sample,Spike_mAbs_inhibitors,3CLpro_inhibitors,RdRp_inhibitors
    ${meta.id},No Mutations,No Mutations,No Mutations
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
