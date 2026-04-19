process CHOPPER {
    tag "$meta.id"
    label 'process_medium'
    //errorStrategy 'ignore'

    container 'quay.io/biocontainers/chopper:0.9.0--hdcf5f25_0'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*fastq") , emit: chopperfastq
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    chopper -q 11  -i $fastq > ${meta.id}_filtered.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: \$(chopper --version 2>&1 | head -n 1)
    END_VERSIONS
    """

    stub:
    """
    cat <<EOF > ${meta.id}_filtered.fastq
    @${meta.id}
    ACGT
    +
    !!!!
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: "stub"
    END_VERSIONS
    """
}
