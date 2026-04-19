# nf-core/sars: Output

## Introduction

This document describes the main outputs produced by the workflow variants in this repository. Exact directory contents depend on the selected `--file` mode.

## Common output areas

### `pipeline_info/`

Run metadata written by Nextflow:

- `execution_report_*.html`
- `execution_timeline_*.html`
- `execution_trace_*.txt`
- `pipeline_dag_*.html`

### `split_fastas/`

Produced by the FASTA-oriented workflows when multi-record FASTA files are split into per-sequence FASTA inputs.

## FASTQ workflow outputs

### Consensus and alignment outputs

The FASTQ workflow generates consensus FASTA files, BAM files, and related ARTIC outputs for each sample.

### `depth/`

Per-position depth summaries derived from the ARTIC alignment workflow, including CSV tables for downstream QC review.

### `primer_metrics/`

Primer database and mismatch summaries used for laboratory QC and reporting. This includes the generated primer database JSON and per-sample mismatch tables.

### `nextclade/`

Nextclade lineage, clade, QC, and mutation summaries converted into report-ready tables.

### Reporting tables

The reporting modules generate CSV summary tables for:

- sequence-level QC metrics
- mutation summaries
- resistance annotation lookups

## FASTA workflow outputs

The FASTA workflow focuses on annotation rather than consensus generation. It produces:

- per-sequence Nextclade outputs
- converted mutation and QC tables
- resistance annotation tables
- aggregated final report tables

## Primercheck workflow outputs

The primercheck workflow produces:

- per-sequence split FASTA files
- generated primer database resources
- primer mismatch tables for each consensus sequence

These outputs are intended for primer performance review and longitudinal assay QC.
