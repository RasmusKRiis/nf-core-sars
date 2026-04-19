process AMPLIGONE {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'

    container 'quay.io/biocontainers/ampligone:1.3.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(fastq)
    path(primerdir)

    output:
    tuple val(meta), path("${meta.id}_primer_cleaned.fastq") , emit: primertrimmedfastq
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    ampligone \
        --input $fastq\
        --output ${meta.id}_primer_cleaned.fastq \
        --reference $primerdir/SARS-CoV-2.reference.fasta \
        --primers $primerdir/SARS-CoV-2.primers.fasta \
        --export-primers ${meta.id}_removed_coordinates.bed \
        --amplicon-type end-to-end

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampligone: $(ampligone --version 2>&1 | head -n 1)
    END_VERSIONS
    """

    stub:
    """
    cat <<EOF > ${meta.id}_primer_cleaned.fastq
    @${meta.id}
    ACGT
    +
    !!!!
    EOF

    touch ${meta.id}_removed_coordinates.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampligone: "stub"
    END_VERSIONS
    """
}
