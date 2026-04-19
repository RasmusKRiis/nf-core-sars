process MAKE_DEPTH_MASK {
    tag { "${meta.id}" }
    publishDir "results/qc", mode: 'copy', overwrite: true
    //errorStrategy 'ignore'

    container 'community.wave.seqera.io/library/bedtools_samtools:2932e857ecf6b5f2'
    input:
    tuple val(meta), path(bam), path(bai)
    val   min_depth

    output:
    tuple val(meta), path("${meta.id}.lowcov.bed"), emit: lowcov
    path("versions.yml"), emit: versions

    script:
    """
    set -euo pipefail

    # 1) Per-base depth
    samtools depth -a ${bam} > ${meta.id}.depth.txt

    # 2) Sites with depth < ${min_depth} -> BED
    awk 'BEGIN{OFS="\\t"} { if (\$3 < ${min_depth}) print \$1, \$2-1, \$2 }' ${meta.id}.depth.txt \
      | bedtools merge > ${meta.id}.lowcov.bed || touch ${meta.id}.lowcov.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: $(samtools --version 2>&1 | head -n 1 | sed 's/^samtools //')
        bedtools: $(bedtools --version 2>&1 | head -n 1 | sed 's/^bedtools v//')
    END_VERSIONS
    """

    stub:
    """
    cat <<EOF > ${meta.id}.lowcov.bed
    MN908947.3	0	10
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "stub"
        bedtools: "stub"
    END_VERSIONS
    """
}
