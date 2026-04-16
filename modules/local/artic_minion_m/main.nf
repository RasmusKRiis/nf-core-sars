process ARTIC_MINION_M {
  tag { "${meta.id}" }
  label 'process_high'
  //errorStrategy 'ignore'
  container 'community.wave.seqera.io/library/artic:1.6.2--d4956cdc155b8612'
  cpus { params.artic_threads }

  input:
    tuple val(meta), path(gp_fastq)
    path  bed
    path  ref

  output:
    tuple val(meta), path("${meta.id}.consensus.fasta"), emit: artic_consensus
    path("${meta.id}.consensus.fasta"), emit: artic_consensus_report
    tuple val(meta), path("${meta.id}.consensus.iupac.fasta"), emit: artic_consensus_iupac
    tuple val(meta), path("${meta.id}.consensus.artic-original.fasta"), emit: artic_consensus_original
    tuple val(meta), path("${meta.id}.consensus.iupac.report.txt"), emit: artic_consensus_iupac_report
    tuple val(meta), path("${meta.id}.normalised.vcf.gz"), path("${meta.id}.normalised.vcf.gz.tbi"), emit: artic_normalised_vcf
    tuple val(meta), path("${meta.id}.normalised.iupac-af.vcf.gz"), path("${meta.id}.normalised.iupac-af.vcf.gz.tbi"), emit: artic_iupac_vcf
    tuple val(meta),
          path("${meta.id}.primertrimmed.rg.sorted.bam"),
          path("${meta.id}.primertrimmed.rg.sorted.bam.bai"),
          emit: artic_bam

  script:
  def MODELOPT = (params.artic_model ? "--model ${params.artic_model}" : "").trim()
  def AMBIGMIN = params.artic_iupac_min_af.toString()
  def AMBIGMAX = params.artic_iupac_max_af.toString()
  def MINDEPTH = params.min_depth.toString()

  def cmd = '''
set -euo pipefail

echo "=== ARTIC_MINION_M: BED/REF mode (consensus-only) ==="
echo "PWD: $(pwd)"
echo "BED: __BED__"
echo "REF: __REF__"
[ -s "__BED__" ] || { echo "ERROR: BED missing/empty: __BED__"; exit 2; }
[ -s "__REF__" ] || { echo "ERROR: REF missing/empty: __REF__"; exit 2; }

# Normalize BED to 7 columns
BED7=bed.normalized.bed
awk -F'\\t' '
  BEGIN { OFS="\\t" }
  /^#/ { print; next }
  NF>=7 { print; next }
  NF==6 { len=$3-$2; if(len<1)len=1; seq=""; for(i=0;i<len;i++)seq=seq "A"; print $1,$2,$3,$4,$5,$6,seq; next }
  NF==5 { s="+"; if (toupper($4) ~ /RIGHT/) s="-"; len=$3-$2; if(len<1)len=1; seq=""; for(i=0;i<len;i++)seq=seq "A"; print $1,$2,$3,$4,$5,s,seq; next }
  { print "ERROR: Unsupported BED line with " NF " columns: " $0 > "/dev/stderr"; exit 11 }
' "__BED__" > "$BED7"

# Clair3 models (best effort)
MODELDIR="$PWD/clair3_models"
mkdir -p "$MODELDIR"
artic_get_models --model-dir "$MODELDIR" || true

# Run ARTIC
artic minion \
  --normalise __NORMALISE__ \
  --threads __THREADS__ \
  --bed "$BED7" \
  --ref "__REF__" \
  --model-dir "$MODELDIR" \
  --min-depth __MINDEPTH__ \
  __MODELOPT__ \
  --read-file __GPFASTQ__ \
  __METAID__

# Find consensus
CONS=""
if [ -f "__METAID__.consensus.fasta" ]; then
  CONS="__METAID__.consensus.fasta"
elif ls __METAID__/*.consensus.fasta >/dev/null 2>&1; then
  for f in __METAID__/*.consensus.fasta; do CONS="$f"; break; done
else
  echo "ERROR: Consensus fasta not found." >&2
  ls -lah || true
  exit 4
fi

# Ensure expected filename
if [ "$CONS" != "__METAID__.consensus.fasta" ]; then
  cp "$CONS" "__METAID__.consensus.fasta"
fi
cp "__METAID__.consensus.fasta" "__METAID__.consensus.artic-original.fasta"

# Locate and normalize the ARTIC normalised VCF as explicit process outputs
NORMVCF=""
if [ -f "__METAID__.normalised.vcf.gz" ]; then
  NORMVCF="__METAID__.normalised.vcf.gz"
elif ls __METAID__/*.normalised.vcf.gz >/dev/null 2>&1; then
  for f in __METAID__/*.normalised.vcf.gz; do NORMVCF="$f"; break; done
fi
if [ -z "$NORMVCF" ]; then
  echo "ERROR: Could not find ARTIC normalised VCF output." >&2
  ls -R || true
  exit 7
fi
if [ "$NORMVCF" != "__METAID__.normalised.vcf.gz" ]; then
  cp "$NORMVCF" "__METAID__.normalised.vcf.gz"
  NORMVCF="__METAID__.normalised.vcf.gz"
fi
if [ -f "${NORMVCF}.tbi" ]; then
  if [ "${NORMVCF}.tbi" != "__METAID__.normalised.vcf.gz.tbi" ]; then
    cp "${NORMVCF}.tbi" "__METAID__.normalised.vcf.gz.tbi"
  fi
else
  tabix -f -p vcf "$NORMVCF"
  if [ "${NORMVCF}.tbi" != "__METAID__.normalised.vcf.gz.tbi" ]; then
    cp "${NORMVCF}.tbi" "__METAID__.normalised.vcf.gz.tbi"
  fi
fi

# Optional IUPAC ambiguity remapping from allele frequencies.
# Keeps ARTIC consensus unchanged and writes iupac to a separate FASTA.
# Uses a BAM pileup branch to capture mixed SNP sites that may be filtered out
# from ARTIC's pass/normalised VCF path.
AMBIG_MIN=__AMBIGMIN__
AMBIG_MAX=__AMBIGMAX__
IUPAC_FASTA="__METAID__.consensus.iupac.fasta"
IUPAC_REPORT="__METAID__.consensus.iupac.report.txt"
{
  echo "sample=__METAID__"
  echo "status=not_run"
  echo "changed_positions=0"
  echo "min_af=$AMBIG_MIN"
  echo "max_af=$AMBIG_MAX"
  echo "min_depth=49"
} > "$IUPAC_REPORT"
# Start IUPAC fasta from ARTIC consensus so downstream consensus is never altered.
cp "__METAID__.consensus.fasta" "$IUPAC_FASTA"
# Default iupac VCF output is a copy of original normalised VCF (overwritten if remapping succeeds)
cp "__METAID__.normalised.vcf.gz" "__METAID__.normalised.iupac-af.vcf.gz"
cp "__METAID__.normalised.vcf.gz.tbi" "__METAID__.normalised.iupac-af.vcf.gz.tbi"
if [ "$(awk "BEGIN{print ($AMBIG_MIN>=0 && $AMBIG_MAX<=1 && $AMBIG_MIN<$AMBIG_MAX)?1:0}")" = "1" ]; then
  IUPAC_BAM=""

  if [ -f "__METAID__.primertrimmed.rg.sorted.bam" ]; then
    IUPAC_BAM="__METAID__.primertrimmed.rg.sorted.bam"
  elif ls __METAID__/*.primertrimmed.rg.sorted.bam >/dev/null 2>&1; then
    for f in __METAID__/*.primertrimmed.rg.sorted.bam; do IUPAC_BAM="$f"; break; done
  fi

  if [ -n "$IUPAC_BAM" ]; then
    echo "Applying IUPAC ambiguity consensus from BAM pileup using AF range [$AMBIG_MIN, $AMBIG_MAX]"
    PILEUP_TXT="__METAID__.pileup.txt"
    # Keep producing a VCF artifact for debugging/output while ambiguity calling is pileup-driven.
    PILEUP_VCF="__METAID__.pileup.calls.vcf.gz"
    bcftools mpileup -f "__REF__" -a FORMAT/AD,FORMAT/DP -Ou "$IUPAC_BAM" | \
      bcftools call -mv -Oz -o "$PILEUP_VCF" || true
    [ -s "$PILEUP_VCF" ] && tabix -f -p vcf "$PILEUP_VCF" || true

    samtools mpileup -aa -A -d 0 -Q 0 -f "__REF__" "$IUPAC_BAM" > "$PILEUP_TXT"
    python3 - "$IUPAC_FASTA" "$PILEUP_TXT" "$AMBIG_MIN" "$AMBIG_MAX" "49" "$IUPAC_FASTA" "$IUPAC_REPORT" <<'PY'
import re
import sys
from collections import Counter

fasta_path, pileup_path, min_af_s, max_af_s, min_depth_s, out_path, report_path = sys.argv[1:]
min_af = float(min_af_s)
max_af = float(max_af_s)
min_depth = int(min_depth_s)

iupac = {
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["G", "C"]): "S",
    frozenset(["A", "T"]): "W",
    frozenset(["G", "T"]): "K",
    frozenset(["A", "C"]): "M",
    frozenset(["C", "G", "T"]): "B",
    frozenset(["A", "G", "T"]): "D",
    frozenset(["A", "C", "T"]): "H",
    frozenset(["A", "C", "G"]): "V",
}

def read_fasta(path):
    header = None
    seq_parts = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is None:
                    header = line
                continue
            seq_parts.append(line.upper())
    if header is None:
        raise RuntimeError(f"No FASTA header in {path}")
    return header, list("".join(seq_parts))

def clean_bases(bases):
    out = []
    i = 0
    n = len(bases)
    while i < n:
        c = bases[i]
        if c == "^":
            i += 2
            continue
        if c == "$":
            i += 1
            continue
        if c in "+-":
            i += 1
            j = i
            while j < n and bases[j].isdigit():
                j += 1
            if j == i:
                continue
            indel_len = int(bases[i:j])
            i = j + indel_len
            continue
        out.append(c)
        i += 1
    return out

header, seq = read_fasta(fasta_path)
changed = 0

with open(pileup_path) as fh:
    for line in fh:
        cols = line.rstrip("\\n").split("\\t")
        if len(cols) < 5:
            continue
        pos = int(cols[1])
        ref = cols[2].upper()
        depth = int(cols[3])
        bases = cols[4]
        if depth < min_depth:
            continue
        counts = Counter()
        for b in clean_bases(bases):
            if b in ".,":  # match reference
                base = ref
            else:
                base = b.upper()
            if base in ("A", "C", "G", "T"):
                counts[base] += 1
        total = sum(counts.values())
        if total < min_depth:
            continue
        selected = sorted(
            [b for b, c in counts.items() if min_af <= (c / total) <= max_af]
        )
        if len(selected) < 2:
            continue
        code = iupac.get(frozenset(selected))
        if not code:
            continue
        idx = pos - 1
        if 0 <= idx < len(seq):
            if seq[idx] != code:
                seq[idx] = code
                changed += 1

with open(out_path, "w") as out:
    out.write(header + "\\n")
    s = "".join(seq)
    for i in range(0, len(s), 80):
        out.write(s[i:i+80] + "\\n")

print(f"IUPAC pileup remap changed_positions={changed}")
with open(report_path, "w") as rep:
    rep.write(f"status=success\\n")
    rep.write(f"changed_positions={changed}\\n")
    rep.write(f"min_af={min_af}\\n")
    rep.write(f"max_af={max_af}\\n")
    rep.write(f"min_depth={min_depth}\\n")
PY
  else
    echo "WARNING: Missing primertrimmed-bam for IUPAC remapping; keeping ARTIC consensus unchanged." >&2
    {
      echo "status=missing_inputs"
      echo "changed_positions=0"
      echo "min_af=$AMBIG_MIN"
      echo "max_af=$AMBIG_MAX"
      echo "min_depth=49"
    } > "$IUPAC_REPORT"
  fi
else
  echo "WARNING: Invalid AF range for IUPAC remapping (min=$AMBIG_MIN max=$AMBIG_MAX); skipping." >&2
  {
    echo "status=invalid_af_range"
    echo "changed_positions=0"
    echo "min_af=$AMBIG_MIN"
    echo "max_af=$AMBIG_MAX"
    echo "min_depth=49"
  } > "$IUPAC_REPORT"
fi

# Make header super simple: >__METAID__
awk -v H=">__METAID__" '/^>/{print H; next} {print}' \
  "__METAID__.consensus.fasta" > "__METAID__.consensus.tmp" && \
  mv "__METAID__.consensus.tmp" "__METAID__.consensus.fasta"
awk -v H=">__METAID__" '/^>/{print H; next} {print}' \
  "__METAID__.consensus.iupac.fasta" > "__METAID__.consensus.iupac.tmp" && \
  mv "__METAID__.consensus.iupac.tmp" "__METAID__.consensus.iupac.fasta"
awk -v H=">__METAID__" '/^>/{print H; next} {print}' \
  "__METAID__.consensus.artic-original.fasta" > "__METAID__.consensus.artic-original.tmp" && \
  mv "__METAID__.consensus.artic-original.tmp" "__METAID__.consensus.artic-original.fasta"

# Grab primer-trimmed BAM for downstream depth analysis
BAM_CANDIDATE=""
if [ -f "__METAID__.primertrimmed.rg.sorted.bam" ]; then
  BAM_CANDIDATE="__METAID__.primertrimmed.rg.sorted.bam"
elif ls __METAID__/*.primertrimmed.rg.sorted.bam >/dev/null 2>&1; then
  for f in __METAID__/*.primertrimmed.rg.sorted.bam; do BAM_CANDIDATE="$f"; break; done
fi

if [ -z "$BAM_CANDIDATE" ]; then
  echo "ERROR: Could not find primertrimmed sorted BAM from ARTIC output." >&2
  ls -R || true
  exit 6
fi

TARGET_BAM="__METAID__.primertrimmed.rg.sorted.bam"
if [ "$BAM_CANDIDATE" != "$TARGET_BAM" ]; then
  cp "$BAM_CANDIDATE" "$TARGET_BAM"
else
  echo "Reusing BAM already at expected path: $TARGET_BAM"
fi

TARGET_BAI="${TARGET_BAM}.bai"
BAM_INDEX_SOURCE="${BAM_CANDIDATE}.bai"
if [ "$BAM_INDEX_SOURCE" != "$TARGET_BAI" ] && [ -f "$BAM_INDEX_SOURCE" ]; then
  cp "$BAM_INDEX_SOURCE" "$TARGET_BAI"
elif [ -f "$TARGET_BAI" ]; then
  echo "Existing BAM index found for $TARGET_BAM"
else
  samtools index -b "$TARGET_BAM"
fi

echo "Final header:"
grep -m1 '^>' "__METAID__.consensus.fasta" || true
'''

  cmd = cmd
    .replace('__BED__', bed.toString())
    .replace('__REF__', ref.toString())
    .replace('__METAID__', meta.id.toString())
    .replace('__GPFASTQ__', gp_fastq.toString())
    .replace('__THREADS__', task.cpus.toString())
    .replace('__NORMALISE__', params.artic_normalise.toString())
    .replace('__AMBIGMIN__', AMBIGMIN)
    .replace('__AMBIGMAX__', AMBIGMAX)
    .replace('__MINDEPTH__', MINDEPTH)
    .replace('__MODELOPT__', MODELOPT)

  return cmd
}
