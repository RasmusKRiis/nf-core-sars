process PRIMER_MASK {
  tag { "${meta.id}" }
  publishDir "results/qc", mode: 'copy', overwrite: true

  container 'community.wave.seqera.io/library/bedtools_samtools:2932e857ecf6b5f2'

  input:
  tuple val(meta), path(lowcov_bed)
  path primer_bed
  val  mask_primer_ends

  output:
  tuple val(meta), path("${meta.id}.mask.bed"), emit: primer_mask
  path("versions.yml"), emit: versions

  when:
  mask_primer_ends

  script:
  """
  set -euo pipefail

  # Clean + merge lowcov + primer coords to strict 3-col BED
  (
    cat ${lowcov_bed} ${primer_bed} 2>/dev/null || cat ${lowcov_bed}
  ) \
  | tr -d '\\r' \
  | awk 'BEGIN{OFS="\\t"} !/^(track|browser)/ && !/^#/ && NF>=3 {print \$1,\$2,\$3}' \
  | sed 's/[[:space:]]\\+\$//' \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  > ${meta.id}.mask.bed

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bedtools: $(bedtools --version 2>&1 | head -n 1 | sed 's/^bedtools v//')
  END_VERSIONS
  """

  stub:
  """
  cat <<EOF > ${meta.id}.mask.bed
  MN908947.3	0	20
  EOF

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bedtools: "stub"
  END_VERSIONS
  """
}
