#!/usr/bin/env python3
"""
Build BioEdit-ready primer alignments against a reference genome.

The script creates two aligned FASTA outputs:
  1) reference-sense: RIGHT primers are reverse-complemented so the aligned
     bases match the reference orientation.
  2) oligo-sense: primer sequences are placed exactly as supplied in the
     primer FASTA, which preserves the physical oligo sequence.

All primer rows are padded with "-" so each sequence spans the full
reference length and can be opened as an alignment in BioEdit.
"""

from __future__ import annotations

import argparse
import pathlib
from dataclasses import dataclass


def load_fasta_sequences(fasta_path: pathlib.Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current_name: str | None = None

    with fasta_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_name = line[1:].split()[0]
                sequences[current_name] = []
                continue
            if current_name is None:
                raise ValueError(f"Malformed FASTA file: {fasta_path}")
            sequences[current_name].append(line.upper())

    return {name: "".join(parts) for name, parts in sequences.items()}


def wrap_fasta_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def reverse_complement(sequence: str) -> str:
    translation = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")
    return sequence.upper().translate(translation)[::-1]


def infer_strand(primer_name: str) -> str:
    upper = primer_name.upper()
    if "_RIGHT" in upper:
        return "-"
    return "+"


@dataclass(frozen=True)
class PrimerRecord:
    chrom: str
    start0: int
    end0: int
    declared_start0: int
    declared_end0: int
    name: str
    pool: str
    strand: str
    oligo_sequence: str
    placement_source: str
    bed_mismatch_count: int

    @property
    def start1(self) -> int:
        return self.start0 + 1

    @property
    def end1(self) -> int:
        return self.end0

    @property
    def reference_sense_sequence(self) -> str:
        if self.strand == "-":
            return reverse_complement(self.oligo_sequence)
        return self.oligo_sequence


def load_primers_from_bed(
    bed_path: pathlib.Path,
    primer_sequences: dict[str, str],
    reference_sequence: str,
) -> list[PrimerRecord]:
    records: list[PrimerRecord] = []

    with bed_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue

            parts = line.split("\t")
            if len(parts) < 4:
                continue

            chrom = parts[0]
            start0 = int(parts[1])
            end0 = int(parts[2])
            name = parts[3]
            pool = parts[4].strip() if len(parts) > 4 else ""
            strand = parts[5].strip() if len(parts) > 5 else infer_strand(name)
            sequence = primer_sequences.get(name)
            if sequence is None:
                raise ValueError(f"Primer {name} exists in BED but not in FASTA.")

            reference_sense_sequence = sequence if strand == "+" else reverse_complement(sequence)
            span_matches = len(sequence) == (end0 - start0)
            bed_mismatch_count = -1
            if span_matches:
                bed_segment = reference_sequence[start0:end0].upper()
                bed_mismatch_count = sum(
                    1 for observed, expected in zip(reference_sense_sequence, bed_segment) if observed != expected
                )
            placement_source = "bed"
            resolved_start0 = start0
            resolved_end0 = end0

            should_relocate = not span_matches
            if span_matches and bed_mismatch_count > 2:
                hit_start0 = reference_sequence.find(reference_sense_sequence)
                should_relocate = hit_start0 != -1
            if should_relocate:
                hit_start0 = reference_sequence.find(reference_sense_sequence)
                if hit_start0 == -1:
                    if not span_matches:
                        raise ValueError(
                            f"Primer {name} has a BED span/sequence length mismatch and "
                            "could not be located in the reference by sequence."
                        )
                else:
                    second_hit = reference_sequence.find(reference_sense_sequence, hit_start0 + 1)
                    if second_hit != -1:
                        raise ValueError(
                            f"Primer {name} requires relocation but maps to multiple "
                            "reference positions, so automatic placement is unsafe."
                        )
                    resolved_start0 = hit_start0
                    resolved_end0 = hit_start0 + len(sequence)
                    placement_source = "reference_match"

            if span_matches and placement_source == "bed" and bed_mismatch_count > 0:
                placement_source = "bed_with_mismatches"

            if not span_matches and placement_source != "reference_match":
                raise ValueError(
                    f"Primer {name} has a BED span/sequence length mismatch and "
                    "could not be safely placed."
                )

            if not span_matches and placement_source == "reference_match":
                bed_mismatch_count = len(sequence)

            if span_matches and bed_mismatch_count > 2 and placement_source != "reference_match":
                raise ValueError(
                    f"Primer {name} has {bed_mismatch_count} mismatches at its BED position "
                    "but could not be safely relocated."
                )

            if placement_source == "reference_match" and span_matches:
                relocated_segment = reference_sequence[resolved_start0:resolved_end0].upper()
                relocated_mismatches = sum(
                    1 for observed, expected in zip(reference_sense_sequence, relocated_segment) if observed != expected
                )
                if relocated_mismatches != 0:
                    raise ValueError(
                        f"Primer {name} was relocated but still does not exactly match the reference."
                    )

            records.append(
                PrimerRecord(
                    chrom=chrom,
                    start0=resolved_start0,
                    end0=resolved_end0,
                    declared_start0=start0,
                    declared_end0=end0,
                    name=name,
                    pool=pool,
                    strand=strand,
                    oligo_sequence=sequence,
                    placement_source=placement_source,
                    bed_mismatch_count=bed_mismatch_count,
                )
            )

    return records


def build_gapped_sequence(reference_length: int, start0: int, sequence: str) -> str:
    if start0 < 0 or (start0 + len(sequence)) > reference_length:
        raise ValueError("Primer coordinates fall outside the reference length.")
    return "-" * start0 + sequence + "-" * (reference_length - start0 - len(sequence))


def write_alignment(
    output_path: pathlib.Path,
    reference_name: str,
    reference_sequence: str,
    records: list[PrimerRecord],
    sequence_getter,
) -> None:
    with output_path.open("w") as handle:
        handle.write(f">{reference_name}\n")
        handle.write(wrap_fasta_sequence(reference_sequence))
        handle.write("\n")

        for record in records:
            placed_sequence = build_gapped_sequence(
                reference_length=len(reference_sequence),
                start0=record.start0,
                sequence=sequence_getter(record),
            )
            header = (
                f"{record.name} "
                f"chrom={record.chrom} start={record.start1} end={record.end1} "
                f"strand={record.strand} pool={record.pool or 'NA'} "
                f"source={record.placement_source}"
            )
            handle.write(f">{header}\n")
            handle.write(wrap_fasta_sequence(placed_sequence))
            handle.write("\n")


def write_manifest(output_path: pathlib.Path, records: list[PrimerRecord]) -> None:
    with output_path.open("w") as handle:
        handle.write(
            "\t".join(
                [
                    "name",
                    "chrom",
                    "declared_start_1based",
                    "declared_end_1based",
                    "start_1based",
                    "end_1based",
                    "strand",
                    "pool",
                    "placement_source",
                    "bed_mismatch_count",
                    "length_nt",
                    "oligo_sequence",
                    "reference_sense_sequence",
                ]
            )
        )
        handle.write("\n")

        for record in records:
            handle.write(
                "\t".join(
                    [
                        record.name,
                        record.chrom,
                        str(record.declared_start0 + 1),
                        str(record.declared_end0),
                        str(record.start1),
                        str(record.end1),
                        record.strand,
                        record.pool or "NA",
                        record.placement_source,
                        str(record.bed_mismatch_count),
                        str(len(record.oligo_sequence)),
                        record.oligo_sequence,
                        record.reference_sense_sequence,
                    ]
                )
            )
            handle.write("\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=pathlib.Path, required=True, help="Reference FASTA.")
    parser.add_argument("--primers", type=pathlib.Path, required=True, help="Primer FASTA.")
    parser.add_argument("--scheme", type=pathlib.Path, required=True, help="Primer BED/scheme.")
    parser.add_argument(
        "--output-prefix",
        type=pathlib.Path,
        required=True,
        help="Prefix for generated files. Example: assets/VMIDT.2.2/SARS-CoV-2.primers.aligned",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    reference_records = load_fasta_sequences(args.reference)
    if len(reference_records) != 1:
        raise ValueError("Reference FASTA must contain exactly one sequence.")

    reference_name, reference_sequence = next(iter(reference_records.items()))
    primer_sequences = load_fasta_sequences(args.primers)
    records = load_primers_from_bed(args.scheme, primer_sequences, reference_sequence)

    write_alignment(
        output_path=pathlib.Path(f"{args.output_prefix}.refsense.fasta"),
        reference_name=reference_name,
        reference_sequence=reference_sequence,
        records=records,
        sequence_getter=lambda record: record.reference_sense_sequence,
    )
    write_alignment(
        output_path=pathlib.Path(f"{args.output_prefix}.oligo.fasta"),
        reference_name=reference_name,
        reference_sequence=reference_sequence,
        records=records,
        sequence_getter=lambda record: record.oligo_sequence,
    )
    write_manifest(
        output_path=pathlib.Path(f"{args.output_prefix}.manifest.tsv"),
        records=records,
    )


if __name__ == "__main__":
    main()
