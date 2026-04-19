process MK_MASK_NO_PRIMER {
  tag { "${meta.id}" }
  publishDir "results/qc", mode: 'copy', overwrite: true

  input:
    tuple val(meta), path(lowcov_bed)
    val  mask_primer_ends

  when:
    !mask_primer_ends

  output:
    tuple val(meta), path("${meta.id}.mask.bed"), emit: no_primer_mask  
    path("versions.yml"), emit: versions

  container 'bash:5.2'
  script:
  """
  set -euo pipefail

  cp ${lowcov_bed} ${meta.id}.mask.bed || touch ${meta.id}.mask.bed

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bash: $(bash --version 2>&1 | head -n 1)
  END_VERSIONS
  """

  stub:
  """
  cat <<EOF > ${meta.id}.mask.bed
  MN908947.3	0	10
  EOF

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bash: "stub"
  END_VERSIONS
  """
}
