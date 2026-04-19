process BCFTOOLS_CONSENSUS {
  tag { "${meta.id}" }
  publishDir "results/consensus", mode: 'copy', overwrite: true
  container 'community.wave.seqera.io/library/bcftools:1.22--a51ee80717c2467e'

  input:
    tuple val(meta), path(vcf), path(mask_bed)
    path  reference

  output:
    tuple val(meta), path("${meta.id}.consensus.fasta"), path("${meta.id}.consensus.report.txt")
    tuple val(meta), path("${meta.id}.consensus.fasta"), emit: bcft_consensus
    path("${meta.id}.consensus.fasta"), emit: bcft_report_consensus
    path("versions.yml"), emit: versions

  script:
  """
  set -euo pipefail

  # Ensure VCF is indexed
  bcftools index -t ${vcf} || true

  # Apply variants; if mask exists/non-empty, apply it in the SAME call
  if [ -s ${mask_bed} ]; then
    bcftools consensus -f ${reference} -m ${mask_bed} ${vcf} > ${meta.id}.consensus.fasta
  else
    bcftools consensus -f ${reference} ${vcf} > ${meta.id}.consensus.fasta
  fi

  printf "sample\\t${meta.id}\\nvcf\\t${vcf}\\nmask\\t${mask_bed}\\n" > ${meta.id}.consensus.report.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bcftools: $(bcftools --version 2>&1 | head -n 1 | sed 's/^bcftools //')
  END_VERSIONS
  """

  stub:
  """
  cat <<EOF > ${meta.id}.consensus.fasta
  >${meta.id}
  ACGTAC
  EOF

  cat <<EOF > ${meta.id}.consensus.report.txt
  sample	${meta.id}
  vcf	stub
  mask	stub
  EOF

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bcftools: "stub"
  END_VERSIONS
  """
}
