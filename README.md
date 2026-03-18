# sarsseq :high_brightness:

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Introduction

This pipeline processes FASTQ files from Nanopore sequencing of SARS-CoV-2, generating consensus sequences and analyzing them for mutations, sequencing statistics, and drug resistance effects. The main steps include:

- Alignment and consensus sequencing with IRMA.
- Consensus sequence analysis with Nextclade.
- Mutation calling
- Generation of a comprehensive report in CSV format.
- Output of sequences in a multiple FASTA file.

The pipeline consist of four different worflows listed bellow:

1) SARS-CoV-2 FASTQ analysis (human)
  Alignment of FASTQ and mutation analysis 
2) SARS-CoV-2 FASTA analysis (human-fasta) (under development)
   Mutation analysis 

## Compatibility

- **Operating System**: Linux
- **Dependencies**: Docker and Nextflow

## Usage

### Sample Sheet Preparation

Prepare a sample sheet (CSV or TSV*) in the `assets` folder with the following format:
* TSV file is not not compulsory

```
PCR-PlatePosition,SequenceID,Barcode,KonsCt
A1*,sampleID,barcodeID,ct-value*
```
*not compulsory

Each row lists a sample to be analyzed. Samples not listed in the sheet will be excluded from the analysis.

### Directory Structure

#### For FASTQ-analysis
Ensure your directory structure is as follows:

```
./
  |-data
         |-barcode3
               |-XXXX_pass_barcode03_XXXX.fastq.gz
               |-YYYY_pass_barcode03_YYYY.fastq.gz
  |-nf-core-sars
               |-assets
                     |-samplesheet.csv
                     |-samplesheet.tsv
               |-...
```

### Running the Pipeline

Navigate to the `nf-core-sars` folder and execute the following command with default parameters:

#### SARS-CoV-2 FASTQ analysis

```bash
nextflow run main.nf -profile docker --runid runid_name  --primerdir primer_folder  --input samplesheet.csv --outdir ../outdir_name
```

### Important Parameters

- `--input` (default: `assets/samplesheet.csv`): Path to the samplesheet.
- `--samplesDir` (default: `../data`): Directory containing the FASTQ files in the structure given above.
- `--primerdir` (default: `assets/V5.4.2/`): Directory containing the primers used during amplification of target region(s).
- `--primer_bed` (optional): BED file describing the primer scheme used for ARTIC, depth analysis and primer QC. If omitted, the pipeline tries to resolve a `.scheme.bed` file from `--primerdir`.
- `--reference` (optional): Reference FASTA for ARTIC consensus. If omitted, the pipeline tries to resolve a `*.reference.fasta` file from `--primerdir`.
- `--primer_fasta` (optional): Explicit path to the FASTA file containing the primer sequences. If omitted, the pipeline tries to resolve `primers.fasta` relative to `--primerdir` or the BED file.
- `--artic_iupac_min_af` (default: `0.30`): Lower AF bound for converting mixed SNP sites to IUPAC ambiguity codes in final ARTIC consensus.
- `--artic_iupac_max_af` (default: `0.70`): Upper AF bound for converting mixed SNP sites to IUPAC ambiguity codes in final ARTIC consensus.
- `--file primercheck-workflow`: Run the FASTA-only primer QC workflow.

All parameters are detailed in the `nextflow.config` file.

## Pipeline Output

The output includes:

- Consensus sequences.
- Depth-per-position CSV tables with amplicon annotations (`results/depth/*_depth_by_position.csv`).
- A primer database JSON plus per-sample primer mismatch matrices suited for reporting (Power BI) in `results/primer_metrics/`.
- Mutation calls.
- Sequencing statistics (coverage, quality parameters).
- Drug resistance effects.
- A report in CSV format.
- A multiple FASTA file of sequences that passed quality filters.

## Primer-only workflow

When you only need primer mismatch statistics for consensus FASTA files, launch the primer-only workflow:

```bash
nextflow run main.nf -profile docker \
    --file primercheck-workflow \
    --fasta /path/to/multi.fasta \
    --primer_bed /path/to/primer.scheme.bed \
    --primer_fasta /path/to/primers.fasta \
    --primer_set_name VM.3 \
    --runid VM3_PRIMER_Q1 \
    --outdir ./primercheck-results
```

This splits each multi-FASTA into per-sequence files, builds the primer database for the provided scheme and writes primer mismatch CSVs to `primer_metrics/` inside `--outdir`.

### SMB helper script

The `wrapper-primercheck.sh` script downloads the requested FASTA/BED/primer FASTA from the N-drive (`//pos1-fhi-svm01.fhi.no/styrt`) via `smbclient` and runs the workflow:

```bash
./wrapper-primercheck.sh \
    -f 2025-02-05_Export.fasta \
    -b VM.3.scheme.bed \
    -P VM.3.primers.fasta \
    -n VM.3 \
    -r VM3_Q1 \
    -o /mnt/tempdata/primercheck-run
```

Store SMB credentials in `~/.smbcreds` (same format as other wrappers).

## Credits

sarsseq was originally written by Rasmus Kopperud Riis.
