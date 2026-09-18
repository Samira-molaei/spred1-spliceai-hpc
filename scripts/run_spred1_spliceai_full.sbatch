#!/bin/bash
#SBATCH --job-name=SPRED1_SpliceAI
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=spred1_spliceai_%j.out
#SBATCH --error=spred1_spliceai_%j.err

set -euo pipefail

WORKDIR="${SPRED1_WORKDIR:-$PWD}"
SIF="${SPLICEAI_SIF:-${WORKDIR}/spliceai-38.sif}"
PY="${PYTHON_BIN:-python3}"
SCRIPT="${BATCH_SCRIPT:-${WORKDIR}/scripts/run_spred1_spliceai_batch.py}"
INPUT="${SPRED1_INPUT:-${WORKDIR}/spred1_hgvs_variants.txt}"
OUTPUT="${SPRED1_OUTPUT:-${WORKDIR}/spred1_spliceai_results_2dp.csv}"

cd "$WORKDIR"

echo "============================================================"
echo "SPRED1 SpliceAI batch"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "============================================================"

"$PY" --version
echo "Input variants: $(grep -cve '^[[:space:]]*$' "$INPUT")"

PORT=$("$PY" - <<'PY'
import socket
s = socket.socket()
s.bind(("", 0))
print(s.getsockname()[1])
s.close()
PY
)

echo "Using local SpliceAI port: $PORT"

export SINGULARITYENV_GENE_SET=basic
export SINGULARITYENV_PORT="$PORT"
export SINGULARITYENV_WORKERS=1

SERVER_LOG="${WORKDIR}/spliceai_server_${SLURM_JOB_ID}.log"

singularity exec --no-home "$SIF" \
  sh -c 'cd / && exec gunicorn \
    --bind ":${PORT}" \
    --workers "${WORKERS}" \
    --threads 1 \
    --timeout 870 \
    --graceful-timeout 30 \
    server:app' \
  > "$SERVER_LOG" 2>&1 &

SERVER_PID=$!
cleanup() { kill "$SERVER_PID" 2>/dev/null || true; }
trap cleanup EXIT

echo "Waiting for SpliceAI server..."
READY=0
for i in $(seq 1 60); do
    if ! kill -0 "$SERVER_PID" 2>/dev/null; then
        echo "ERROR: SpliceAI server exited during startup."
        tail -100 "$SERVER_LOG"
        exit 1
    fi

    if curl -fsS \
      "http://localhost:${PORT}/spliceai/?hg=38&variant=15-38253191-C-G&bc=basic&distance=500&mask=0" \
      >/dev/null 2>&1; then
        READY=1
        break
    fi
    sleep 5
done

if [ "$READY" -ne 1 ]; then
    echo "ERROR: SpliceAI server did not become ready."
    tail -100 "$SERVER_LOG"
    exit 1
fi

export SPLICEAI_API="http://localhost:${PORT}/spliceai/"

echo "SpliceAI server is ready."
echo "API: $SPLICEAI_API"
echo "Settings: hg38 / GENCODE basic / distance=500 / mask=0"
echo "Transcript: ENST00000299084.9 (SPRED1 MANE Select)"
echo "Output display: 2 decimal places"

"$PY" -u "$SCRIPT" "$INPUT" "$OUTPUT" --resume

echo "============================================================"
echo "FINAL QC"
echo "============================================================"
ls -lh "$OUTPUT"
wc -l "$OUTPUT"

"$PY" - "$OUTPUT" <<'PY'
import csv, sys
from collections import Counter
with open(sys.argv[1], newline="") as f:
    c = Counter(row["status"] for row in csv.DictReader(f))
print("Status counts:", dict(c))
PY
