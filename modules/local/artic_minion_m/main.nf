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
    tuple val(meta),
          path("${meta.id}.primertrimmed.rg.sorted.bam"),
          path("${meta.id}.primertrimmed.rg.sorted.bam.bai"),
          emit: artic_bam

  script:
  def MODELOPT = (params.artic_model ? "--model ${params.artic_model}" : "").trim()
  def AMBIGMIN = params.artic_iupac_min_af.toString()
  def AMBIGMAX = params.artic_iupac_max_af.toString()

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
  --min-depth 14 \
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

# Optional IUPAC ambiguity remapping from allele frequencies.
# Keeps original ARTIC consensus as backup and rewrites consensus.fasta in-place.
AMBIG_MIN=__AMBIGMIN__
AMBIG_MAX=__AMBIGMAX__
if [ "$(awk "BEGIN{print ($AMBIG_MIN>=0 && $AMBIG_MAX<=1 && $AMBIG_MIN<$AMBIG_MAX)?1:0}")" = "1" ]; then
  PRECONS=""
  NORMVCF=""
  MASK=""

  if [ -f "__METAID__.preconsensus.fasta" ]; then
    PRECONS="__METAID__.preconsensus.fasta"
  elif ls __METAID__/*.preconsensus.fasta >/dev/null 2>&1; then
    for f in __METAID__/*.preconsensus.fasta; do PRECONS="$f"; break; done
  fi

  if [ -f "__METAID__.normalised.vcf.gz" ]; then
    NORMVCF="__METAID__.normalised.vcf.gz"
  elif ls __METAID__/*.normalised.vcf.gz >/dev/null 2>&1; then
    for f in __METAID__/*.normalised.vcf.gz; do NORMVCF="$f"; break; done
  fi

  if [ -f "__METAID__.coverage_mask.txt" ]; then
    MASK="__METAID__.coverage_mask.txt"
  elif ls __METAID__/*.coverage_mask.txt >/dev/null 2>&1; then
    for f in __METAID__/*.coverage_mask.txt; do MASK="$f"; break; done
  fi

  if [ -n "$PRECONS" ] && [ -n "$NORMVCF" ] && [ -n "$MASK" ]; then
    echo "Applying IUPAC ambiguity consensus using AF range [$AMBIG_MIN, $AMBIG_MAX]"
    cp "__METAID__.consensus.fasta" "__METAID__.consensus.artic-original.fasta"

    AF_EXPR=""
    if bcftools view -h "$NORMVCF" | grep -q '^##FORMAT=<ID=AF,'; then
      AF_EXPR="FMT/AF>=$AMBIG_MIN && FMT/AF<=$AMBIG_MAX"
    elif bcftools view -h "$NORMVCF" | grep -q '^##INFO=<ID=AF,'; then
      AF_EXPR="INFO/AF>=$AMBIG_MIN && INFO/AF<=$AMBIG_MAX"
    elif bcftools view -h "$NORMVCF" | grep -q '^##FORMAT=<ID=AD,'; then
      AF_EXPR="FMT/AD[1]/(FMT/AD[0]+FMT/AD[1])>=$AMBIG_MIN && FMT/AD[1]/(FMT/AD[0]+FMT/AD[1])<=$AMBIG_MAX"
    fi

    if [ -n "$AF_EXPR" ]; then
      AMBIG_VCF="__METAID__.normalised.iupac-af.vcf.gz"
      bcftools +setGT "$NORMVCF" -Oz -o "$AMBIG_VCF" -- -t q -n c:'0/1' -i "$AF_EXPR"
      tabix -f -p vcf "$AMBIG_VCF"
      bcftools consensus --iupac-codes -f "$PRECONS" "$AMBIG_VCF" -m "$MASK" -o "__METAID__.consensus.fasta"
    else
      echo "WARNING: Could not derive AF expression (no AF/AD fields). Keeping ARTIC consensus unchanged." >&2
    fi
  else
    echo "WARNING: Missing preconsensus/normalised-vcf/mask for IUPAC remapping; keeping ARTIC consensus unchanged." >&2
  fi
else
  echo "WARNING: Invalid AF range for IUPAC remapping (min=$AMBIG_MIN max=$AMBIG_MAX); skipping." >&2
fi

# Make header super simple: >__METAID__
awk -v H=">__METAID__" '/^>/{print H; next} {print}' \
  "__METAID__.consensus.fasta" > "__METAID__.consensus.tmp" && \
  mv "__METAID__.consensus.tmp" "__METAID__.consensus.fasta"

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
    .replace('__MODELOPT__', MODELOPT)

  return cmd
}
