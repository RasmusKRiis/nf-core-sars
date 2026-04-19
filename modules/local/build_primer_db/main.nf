process BUILD_PRIMER_DB {
    tag "primer_db"
    label 'process_single'
    publishDir "${params.outdir ?: 'results'}/primer_metrics", mode: 'copy', overwrite: true

    container 'docker.io/rasmuskriis/nextclade-python'
    containerOptions = "-v ${baseDir}/bin:/project-bin"

    input:
    path primer_bed
    path primer_fasta
    val primer_set_name

    output:
    path("${primer_set_name}.primer_db.json"), emit: primer_db
    path("versions.yml"), emit: versions

    script:
    """
    set -euo pipefail

    python3 /project-bin/primer_metrics.py build-db \
        --primer-bed ${primer_bed} \
        --primer-fasta ${primer_fasta} \
        --primer-set-name ${primer_set_name} \
        --output ${primer_set_name}.primer_db.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cat <<EOF > ${primer_set_name}.primer_db.json
    {
      "primer_set": "${primer_set_name}",
      "primers": {},
      "amplicons": {}
    }
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
