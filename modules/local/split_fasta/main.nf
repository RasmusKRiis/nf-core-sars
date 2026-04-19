process SPLIT_FASTA {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/rasmuskriis/nextclade-python'
    containerOptions = "-v ${baseDir}/bin:/project-bin"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('split/*.fasta'), emit: split_fastas
    path('versions.yml'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    mkdir -p split

    python3 - "${fasta}" <<'PY'
import re
import sys
from pathlib import Path

fasta_path = Path(sys.argv[1])
outdir = Path("split")
outdir.mkdir(exist_ok=True)

records = []
header = None
sequence_lines = []

with fasta_path.open() as handle:
    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence_lines)))
            header = line[1:].strip()
            sequence_lines = []
        else:
            sequence_lines.append(line)

if header is not None:
    records.append((header, "".join(sequence_lines)))

if not records:
    raise SystemExit(f"No FASTA records found in {fasta_path}")

def sanitize(header_text):
    base = header_text.split()[0].strip()
    base = re.sub(r"[^A-Za-z0-9._-]+", "_", base)
    base = base.strip("._-")
    return base or "sequence"

counts = {}

for index, (record_header, sequence) in enumerate(records, start=1):
    base_name = sanitize(record_header)
    counts[base_name] = counts.get(base_name, 0) + 1
    if counts[base_name] > 1:
        base_name = f"{base_name}_{counts[base_name]}"

    output_path = outdir / f"{base_name}.fasta"
    with output_path.open("w") as handle:
        handle.write(f">{record_header}\\n")
        for offset in range(0, len(sequence), 80):
            handle.write(sequence[offset:offset + 80] + "\\n")
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
