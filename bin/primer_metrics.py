#!/usr/bin/env python3
"""
Primer utilities for nf-core/sarscovseq:

 1) build-db: generate a JSON primer database mirroring the primer set used.
 2) mismatch: compute primer mismatch metrics for each sample/primer pair.
"""

import argparse
import csv
import json
import pathlib
import re
from typing import Dict, List, Tuple


AMBIGUITY_CODES = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "U": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}


def load_fasta_sequences(fasta_path: pathlib.Path) -> Dict[str, str]:
    seqs = {}
    current = None
    with fasta_path.open("r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            else:
                if current is None:
                    raise ValueError(f"FASTA file {fasta_path} is malformed.")
                seqs[current].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def derive_amplicon_name(primer_name: str) -> str:
    if not primer_name:
        return "UNKNOWN"
    match = re.split(r"_(LEFT|RIGHT)", primer_name, flags=re.IGNORECASE)
    if match:
        return match[0]
    return primer_name


def primer_name_aliases(primer_name: str) -> List[str]:
    """
    Generate alternative FASTA header candidates for a BED primer name.
    Handles common naming differences, for example:
      SARS-CoV-2_10_LEFT_2 <-> SARSCoV_amplicon_10_LEFT_alt1
    """
    aliases: List[str] = [primer_name]

    normalized = primer_name.replace("SARS-CoV-2", "SARSCoV")
    if normalized != primer_name:
        aliases.append(normalized)

    m = re.match(r"^SARS-CoV-2_(\d+)_(LEFT|RIGHT)_(\d+)$", primer_name)
    if m:
        amplicon, side, variant = m.groups()
        base = f"SARSCoV_amplicon_{amplicon}_{side}"
        aliases.append(base)
        try:
            variant_i = int(variant)
            if variant_i > 1:
                aliases.append(f"{base}_alt{variant_i - 1}")
        except ValueError:
            pass

    # De-duplicate while preserving order.
    deduped = []
    seen = set()
    for alias in aliases:
        if alias not in seen:
            seen.add(alias)
            deduped.append(alias)
    return deduped


def build_primer_db(primer_bed: pathlib.Path, primer_fasta: pathlib.Path, primer_set_name: str, output: pathlib.Path) -> None:
    primer_sequences = load_fasta_sequences(primer_fasta)
    primer_sequences_ci = {k.upper(): v for k, v in primer_sequences.items()}
    primers: Dict[str, Dict] = {}
    amplicons: Dict[str, Dict] = {}

    with primer_bed.open("r") as bed:
        for raw_line in bed:
            line = raw_line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            name = parts[3]
            pool = parts[4].strip() if len(parts) > 4 else ""
            strand = parts[5].strip() if len(parts) > 5 else ""
            if strand not in {"+", "-"}:
                upper_name = name.upper()
                if "_RIGHT" in upper_name:
                    strand = "-"
                else:
                    strand = "+"
            seq_from_bed = parts[6].strip().upper() if len(parts) > 6 else ""
            sequence = seq_from_bed
            if not sequence:
                for alias in primer_name_aliases(name):
                    sequence = primer_sequences.get(alias, "")
                    if sequence:
                        break
                    sequence = primer_sequences_ci.get(alias.upper(), "")
                    if sequence:
                        break
            if not sequence:
                raise ValueError(f"Missing primer sequence for {name}. Provide a FASTA containing all primers.")

            amplicon_id = derive_amplicon_name(name)
            start_1 = start + 1
            end_1 = max(start_1, end)
            primers[name] = {
                "name": name,
                "sequence": sequence,
                "chrom": chrom,
                "start": start_1,
                "end": end_1,
                "strand": strand or "+",
                "pool": pool,
                "amplicon": amplicon_id,
                "length": len(sequence),
            }

            amp_key = f"{chrom}:{amplicon_id}"
            entry = amplicons.setdefault(
                amp_key,
                {
                    "amplicon": amplicon_id,
                    "chrom": chrom,
                    "start": start_1,
                    "end": end_1,
                    "pools": set(),
                },
            )
            entry["start"] = min(entry["start"], start_1)
            entry["end"] = max(entry["end"], end_1)
            if pool:
                entry["pools"].add(pool)

    formatted_amplicons = {}
    for amp_key, entry in amplicons.items():
        formatted_amplicons[entry["amplicon"]] = {
            "chrom": entry["chrom"],
            "start": entry["start"],
            "end": entry["end"],
            "pools": sorted(entry["pools"]),
            "length": entry["end"] - entry["start"] + 1,
        }

    db = {
        "primer_set": primer_set_name,
        "primer_bed": str(primer_bed),
        "primer_fasta": str(primer_fasta),
        "primers": primers,
        "amplicons": formatted_amplicons,
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        json.dump(db, handle, indent=2)


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")
    return seq.upper().translate(comp)[::-1]


def bases_match(primer_base: str, template_base: str) -> bool:
    """
    Compare primer and template bases with IUPAC support.
    Unknown template bases should be penalized: treating template N as a match can
    produce false "perfect" hits in low-quality/ambiguous regions.
    """
    t_base = template_base.upper()
    if t_base == "N" or t_base == "-":
        return False
    p_set = AMBIGUITY_CODES.get(primer_base.upper(), {primer_base.upper()})
    t_set = AMBIGUITY_CODES.get(t_base, {t_base})
    return bool(p_set & t_set)


def mismatch_metrics(primer_seq: str, template_seq: str) -> Tuple[int, List[int]]:
    mismatches = []
    for idx, (a, b) in enumerate(zip(primer_seq, template_seq), start=1):
        if not bases_match(a, b):
            mismatches.append(idx)
    if len(template_seq) < len(primer_seq):
        span = range(len(template_seq) + 1, len(primer_seq) + 1)
        mismatches.extend(span)
    return len(mismatches), mismatches


def load_primer_db(path: pathlib.Path) -> Dict:
    with path.open("r") as handle:
        return json.load(handle)


def load_consensus_sequence(path: pathlib.Path) -> str:
    seqs = load_fasta_sequences(path)
    if not seqs:
        raise ValueError(f"No sequences found in {path}")
    # use the first record
    return next(iter(seqs.values()))


def find_best_alignment(primer_seq: str, window_seq: str) -> Tuple[int, List[int], int, str]:
    """Slide the primer across the window and return the alignment with the fewest mismatches."""
    primer_len = len(primer_seq)
    if primer_len == 0:
        return 0, [], 0, ""
    if len(window_seq) < primer_len:
        window_seq = window_seq + ("N" * (primer_len - len(window_seq)))
    max_offset = len(window_seq) - primer_len
    best = (primer_len + 1, [], 0, window_seq[:primer_len])
    for offset in range(max_offset + 1):
        candidate = window_seq[offset : offset + primer_len]
        mismatch_count, mismatch_positions = mismatch_metrics(primer_seq, candidate)
        if mismatch_count < best[0]:
            best = (mismatch_count, mismatch_positions, offset, candidate)
            if mismatch_count == 0:
                break
    return best


def alignment_coordinates(
    offset: int,
    primer_len: int,
    strand: str,
    window_start: int,
    window_end: int,
) -> Tuple[int, int]:
    """
    Convert an alignment offset (in the oriented search sequence) to 1-based
    consensus coordinates on the original strand.
    """
    if strand == "-":
        window_len = max(0, window_end - window_start)
        start0 = window_start + (window_len - (offset + primer_len))
        end0 = start0 + primer_len
    else:
        start0 = window_start + offset
        end0 = start0 + primer_len
    return start0 + 1, end0


def write_mismatch_matrix(
    consensus_fasta: pathlib.Path,
    primer_db_path: pathlib.Path,
    sample_id: str,
    output: pathlib.Path,
    flank: int = 30,
    run_id: str = "Unknown",
    global_fallback_mismatches: int = 4,
    edge_buffer: int = 3,
    global_min_distance_gain: int = 10,
) -> None:
    primer_db = load_primer_db(primer_db_path)
    consensus = load_consensus_sequence(consensus_fasta)
    primer_set_name = primer_db.get("primer_set", "")
    primers = primer_db.get("primers", {})

    rows = []
    for primer_name, entry in primers.items():
        primer_seq = entry["sequence"].upper()
        start = int(entry["start"])
        end = int(entry["end"])
        flank_len = max(0, flank)
        window_start = max(0, start - 1 - flank_len)
        window_end = min(len(consensus), end + flank_len)
        template_window = consensus[window_start:window_end]
        if entry.get("strand", "+") == "-":
            oriented_window = reverse_complement(template_window)
        else:
            oriented_window = template_window

        local_mismatch_count, local_mismatch_positions, local_best_offset, local_best_segment = find_best_alignment(
            primer_seq, oriented_window
        )
        local_max_offset = max(0, len(oriented_window) - len(primer_seq))
        near_left_edge = local_best_offset <= edge_buffer
        near_right_edge = local_best_offset >= max(0, local_max_offset - edge_buffer)
        local_is_edge_hit = near_left_edge or near_right_edge
        needs_global_fallback = (
            local_mismatch_count >= max(0, global_fallback_mismatches)
            or (local_is_edge_hit and local_mismatch_count > 0)
        )

        best_mismatch_count = local_mismatch_count
        best_mismatch_positions = local_mismatch_positions
        best_offset = local_best_offset
        best_segment = local_best_segment
        best_window_start = window_start
        best_window_end = window_end
        search_mode = "local"

        if needs_global_fallback:
            if entry.get("strand", "+") == "-":
                oriented_global = reverse_complement(consensus)
            else:
                oriented_global = consensus
            (
                global_mismatch_count,
                global_mismatch_positions,
                global_best_offset,
                global_best_segment,
            ) = find_best_alignment(primer_seq, oriented_global)

            local_loc_start, _ = alignment_coordinates(
                local_best_offset, len(primer_seq), entry.get("strand", "+"), window_start, window_end
            )
            global_loc_start, _ = alignment_coordinates(
                global_best_offset, len(primer_seq), entry.get("strand", "+"), 0, len(consensus)
            )
            local_distance = abs(local_loc_start - start)
            global_distance = abs(global_loc_start - start)

            choose_global = (
                (global_mismatch_count < local_mismatch_count)
                or (
                    global_mismatch_count == local_mismatch_count
                    and (global_distance + max(0, global_min_distance_gain)) <= local_distance
                )
            )

            if choose_global:
                best_mismatch_count = global_mismatch_count
                best_mismatch_positions = global_mismatch_positions
                best_offset = global_best_offset
                best_segment = global_best_segment
                best_window_start = 0
                best_window_end = len(consensus)
                search_mode = "global_fallback"

        located_start, located_end = alignment_coordinates(
            best_offset, len(primer_seq), entry.get("strand", "+"), best_window_start, best_window_end
        )
        distance_from_expected = abs(located_start - start)
        percent_identity = (
            ((len(primer_seq) - best_mismatch_count) / len(primer_seq)) * 100 if primer_seq else 0.0
        )

        rows.append(
            {
                "Run_ID": run_id,
                "Sample_ID": sample_id,
                "Primer_Set": primer_set_name,
                "Primer_Name": primer_name,
                "Amplicon": entry.get("amplicon", ""),
                "Pool": entry.get("pool", ""),
                "Chromosome": entry.get("chrom", ""),
                "Start": start,
                "End": end,
                "Strand": entry.get("strand", "+"),
                "Primer_Length": len(primer_seq),
                "Mismatches": best_mismatch_count,
                "Percent_Identity": f"{percent_identity:.2f}",
                "Mismatch_Positions": ";".join(map(str, best_mismatch_positions)),
                "Primer_Sequence": primer_seq,
                "Template_Sequence": best_segment,
                "Alignment_Offset": best_offset,
                "Located_Start": located_start,
                "Located_End": located_end,
                "Distance_From_Expected": distance_from_expected,
                "Search_Mode": search_mode,
            }
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "Run_ID",
        "Sample_ID",
        "Primer_Set",
        "Primer_Name",
        "Amplicon",
        "Pool",
        "Chromosome",
        "Start",
        "End",
        "Strand",
        "Primer_Length",
        "Mismatches",
        "Percent_Identity",
        "Mismatch_Positions",
        "Primer_Sequence",
        "Template_Sequence",
        "Alignment_Offset",
        "Located_Start",
        "Located_End",
        "Distance_From_Expected",
        "Search_Mode",
    ]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Primer helper utilities.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser("build-db", help="Create a primer JSON database.")
    build.add_argument("--primer-bed", required=True, help="Primer BED file.")
    build.add_argument("--primer-fasta", required=True, help="Primer FASTA file.")
    build.add_argument("--primer-set-name", required=True, help="Name for this primer set.")
    build.add_argument("--output", required=True, help="Output JSON path.")

    mismatch = subparsers.add_parser("mismatch", help="Generate primer mismatch matrix for a sample.")
    mismatch.add_argument("--consensus", required=True, help="Consensus FASTA for the sample.")
    mismatch.add_argument("--primer-db", required=True, help="Primer JSON database.")
    mismatch.add_argument("--sample-id", required=True, help="Sample identifier.")
    mismatch.add_argument("--output", required=True, help="Output CSV path.")
    mismatch.add_argument(
        "--flank",
        type=int,
        default=80,
        help="Number of bases to extend on each side of the primer coordinates when extracting the consensus window.",
    )
    mismatch.add_argument(
        "--run-id",
        type=str,
        default="Unknown",
        help="Run identifier to embed into the mismatch output.",
    )
    mismatch.add_argument(
        "--global-fallback-mismatches",
        type=int,
        default=4,
        help="Trigger global search when local mismatch count is greater than or equal to this threshold.",
    )
    mismatch.add_argument(
        "--edge-buffer",
        type=int,
        default=3,
        help="Trigger global search when local best hit is close to the window edge.",
    )
    mismatch.add_argument(
        "--global-min-distance-gain",
        type=int,
        default=10,
        help="For tied mismatch counts, require this much distance improvement before switching to global.",
    )

    return parser.parse_args()


def main():
    args = parse_args()
    if args.command == "build-db":
        build_primer_db(
            pathlib.Path(args.primer_bed),
            pathlib.Path(args.primer_fasta),
            args.primer_set_name,
            pathlib.Path(args.output),
        )
        print(f"[primer-db] Wrote primer database to {args.output}")
    elif args.command == "mismatch":
        write_mismatch_matrix(
            pathlib.Path(args.consensus),
            pathlib.Path(args.primer_db),
            args.sample_id,
            pathlib.Path(args.output),
            flank=args.flank,
            run_id=args.run_id,
            global_fallback_mismatches=args.global_fallback_mismatches,
            edge_buffer=args.edge_buffer,
            global_min_distance_gain=args.global_min_distance_gain,
        )
        print(f"[primer-mismatch] Wrote mismatch matrix to {args.output}")


if __name__ == "__main__":
    main()
