# nf-core/sars: Usage

## Introduction

This repository currently exposes three workflow modes through the `--file` parameter:

- `fastq-workflow`: end-to-end Nanopore FASTQ processing with consensus generation, depth analysis, primer QC, Nextclade summaries, and reporting.
- `fasta-workflow`: FASTA-only analysis for Nextclade annotation and resistance reporting.
- `primercheck-workflow`: primer mismatch QC for consensus FASTA files.

## FASTQ workflow samplesheet

The FASTQ workflow expects a CSV samplesheet supplied with `--input`. The required columns are:

| Column | Description |
| --- | --- |
| `SequenceID` | Sample identifier used in reports and downstream outputs. |
| `Barcode` | Barcode directory name used to find demultiplexed FASTQ files under `--samplesDir`. |
| `KonsCt` | Optional Ct value retained in the metadata sheet. |

Example:

```csv
SequenceID,Barcode,KonsCt
sample_001,barcode01,23.1
sample_002,barcode02,27.8
```

`--samplesDir` should contain demultiplexed barcode folders like:

```text
data/
├── barcode01/
│   ├── runA_pass_barcode01_000.fastq.gz
│   └── runA_pass_barcode01_001.fastq.gz
└── barcode02/
    └── runA_pass_barcode02_000.fastq.gz
```

## Running the FASTQ workflow

```bash
nextflow run . \
    -profile docker \
    --file fastq-workflow \
    --input assets/samplesheet.csv \
    --samplesDir ../data \
    --primerdir assets/primer_schemes/ncov-2019_midnight/v3.0.0 \
    --outdir results
```

Key parameters:

- `--primerdir`: directory containing the primer scheme resources used in the lab.
- `--primer_bed`: explicit BED file override for the primer scheme.
- `--reference`: explicit reference FASTA override.
- `--primer_fasta`: explicit primer FASTA override.
- `--runid`: run identifier shown in the generated reports.
- `--seq_instrument`: sequencing instrument label shown in reports.

## Running the FASTA workflow

Use the FASTA workflow when consensus sequences already exist and only downstream annotation and reporting are needed.

```bash
nextflow run . \
    -profile docker \
    --file fasta-workflow \
    --fasta "/path/to/*.fasta" \
    --runid RUN_001 \
    --outdir fasta-results
```

## Running the primercheck workflow

The primer-check workflow builds the primer database for a scheme and produces mismatch metrics for every sequence in the supplied FASTA file(s).

```bash
nextflow run . \
    -profile docker \
    --file primercheck-workflow \
    --fasta /path/to/consensus.fasta \
    --primer_bed /path/to/primer.scheme.bed \
    --primer_fasta /path/to/primers.fasta \
    --primer_set_name VM.3 \
    --runid VM3_Q1 \
    --outdir primercheck-results
```

## Profiles

The pipeline ships with the standard nf-core execution profiles, including `docker`, `singularity`, `apptainer`, `podman`, `conda`, `mamba`, `test`, and `test_full`.

## Running with parameter files

You can place pipeline parameters in a YAML or JSON file and pass them with `-params-file`:

```bash
nextflow run . -profile docker -params-file params.yaml
```

This is the preferred way to store reusable run settings. Avoid using `-c` for pipeline parameters; use it only for infrastructure-level Nextflow configuration.
