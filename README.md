<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-sars_logo_dark.png">
    <img alt="nf-core/sars" src="docs/images/nf-core-sars_logo_light.png">
  </picture>
</h1>

[![GitHub Actions nf-test Status](https://github.com/RasmusKRiis/nf-core-sars/actions/workflows/nf-test.yml/badge.svg)](https://github.com/RasmusKRiis/nf-core-sars/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/RasmusKRiis/nf-core-sars/actions/workflows/linting.yml/badge.svg)](https://github.com/RasmusKRiis/nf-core-sars/actions/workflows/linting.yml)
[![Nextflow](https://img.shields.io/badge/version-%E2%89%A523.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.2)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

## Introduction

**nf-core/sars** is a SARS-CoV-2 analysis pipeline for Nanopore amplicon sequencing and downstream consensus review workflows. The repository currently contains three entry points:

1. `fastq-workflow` for demultiplexed Nanopore FASTQ inputs, consensus generation, depth analysis, primer mismatch metrics, Nextclade summaries, and resistance annotation.
2. `fasta-workflow` for FASTA-only Nextclade analysis and report generation.
3. `primercheck-workflow` for primer mismatch QC on consensus FASTA files.

The pipeline produces consensus FASTA files, depth tables, primer QC outputs, mutation summaries, resistance annotations, and run-level reports.

## Usage

The full usage guide is in [docs/usage.md](docs/usage.md). The minimal FASTQ workflow command is:

```bash
nextflow run . \
    -profile docker \
    --file fastq-workflow \
    --input assets/samplesheet.csv \
    --samplesDir ../data \
    --primerdir assets/primer_schemes/ncov-2019_midnight/v3.0.0 \
    --outdir results
```

The FASTA-only and primer-check workflows are documented in [docs/usage.md](docs/usage.md) as well.

## Pipeline Output

The output guide is in [docs/output.md](docs/output.md). Typical result areas include:

- `pipeline_info/` for Nextflow execution metadata.
- `depth/` for per-position coverage summaries.
- `primer_metrics/` for primer mismatch metrics and generated primer databases.
- `nextclade/` and report tables for lineage, mutation, and resistance summaries.

## Credits

nf-core/sars was originally written by Rasmus Kopperud Riis.

## Contributions and Support

Please see [`.github/CONTRIBUTING.md`](.github/CONTRIBUTING.md) for contribution guidance. For nf-core specific discussion, the canonical community channel is the [nf-core Slack `#sars` channel](https://nfcore.slack.com/channels/sars).

## Citations

Tool and framework citations are listed in [CITATIONS.md](CITATIONS.md).
