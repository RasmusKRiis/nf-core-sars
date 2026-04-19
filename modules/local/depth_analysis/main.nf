
process DEPTH_ANALYSIS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir ?: 'results'}/depth", mode: 'copy', overwrite: true
    errorStrategy 'ignore'

    container 'docker.io/rasmuskriis/nextclade-python'
    containerOptions = "-v ${baseDir}/bin:/project-bin"

    input:
    tuple val(meta), path(bam), path(bai)
    path primer_bed

    output:
    tuple val(meta), path("${meta.id}_depth_by_position.csv"), emit: depth_report
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    python3 /project-bin/depth_analysis.py \
        --bam ${bam} \
        --sample-id ${meta.id} \
        --primer-bed ${primer_bed} \
        --output ${meta.id}_depth_by_position.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cat <<EOF > ${meta.id}_depth_by_position.csv
    Sample_ID,Position,Depth
    ${meta.id},1,100
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
