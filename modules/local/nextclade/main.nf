process NEXTCLADE {
    tag "${meta.id}"
    label 'process_single'
    errorStrategy 'ignore'

    container 'docker.io/nextstrain/nextclade:latest'
    containerOptions = "-v ${baseDir}/bin:/project-bin"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_nextclade.csv"), emit: nextclade_csv, optional: true
    path('versions.yml'), emit: versions, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    dataset_dir="${meta.id}_nextclade_dataset"
    output_dir="${meta.id}_nextclade_output"

    nextclade dataset get \
        --name 'nextstrain/sars-cov-2/wuhan-hu-1/orfs' \
        --output-dir "\${dataset_dir}"

    nextclade run \
        --input-dataset "\${dataset_dir}" \
        --output-all="\${output_dir}" \
        ${fasta}

    csv_file=\$(find "\${output_dir}" -maxdepth 1 -type f -name '*.csv' | head -n 1)
    if [[ -z "\${csv_file}" ]]; then
        echo "ERROR: Nextclade did not produce a CSV output for ${meta.id}" >&2
        exit 1
    fi

    mv "\${csv_file}" "${meta.id}_nextclade.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(nextclade --version 2>&1 | sed -E 's/^[^0-9]*([0-9].*)$/\\1/')
    END_VERSIONS
    """
}
