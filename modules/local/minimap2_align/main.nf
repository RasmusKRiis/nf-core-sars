process MINIMAP2_ALIGN {
  tag { "${meta.id}" }
  publishDir "results/aln", mode: 'copy', overwrite: true

  container 'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c'
  cpus 4
  memory '4 GB'

  input:
  tuple val(meta), path(reads)    // accepts single FASTQ or list
  path  reference

  output:
  tuple val(meta),
        path("${meta.id}.bam"),
        path("${meta.id}.bam.bai"),
        emit: minimap2
  path("versions.yml"), emit: versions

  script:
  def readList = reads instanceof List ? reads.collect{ it.toString() } : [reads.toString()]
  """
  set -euo pipefail

  minimap2 -ax map-ont -t ${task.cpus} ${reference} ${readList.join(' ')} \
    | samtools sort -@ ${task.cpus} -o ${meta.id}.bam
  # force BAI index (SARS-CoV-2 is tiny, so BAI is fine)
  samtools index -b ${meta.id}.bam

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      minimap2: $(minimap2 --version 2>&1 | head -n 1)
      samtools: $(samtools --version 2>&1 | head -n 1 | sed 's/^samtools //')
  END_VERSIONS
  """

  stub:
  """
  touch ${meta.id}.bam
  touch ${meta.id}.bam.bai

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      minimap2: "stub"
      samtools: "stub"
  END_VERSIONS
  """
}
