#!/usr/bin/env bash
set -euo pipefail

# Batch runner launcher for Taguchi L9 inputs.
# Edit values below, save, then run: ./run_taguchi_l9.sh

PYTHON_BIN="${PYTHON_BIN:-python}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_CSV="${SCRIPT_DIR}/data/taguchi_l9_inputs.csv"
OUTPUT_DIR="${SCRIPT_DIR}/output/taguchi_l9"
OUTPUT_CSV="${OUTPUT_DIR}/taguchi_l9_predictions.csv"
MODEL_JSON=""

# Observation settings used for percentile range output.
SEM_COUNT=10
AGING_HOURS=0.0
AGGLOMERATION_RATE=0.08
GEO_SIGMA=1.35

CMD=(
  "${PYTHON_BIN}" "${SCRIPT_DIR}/run_taguchi_l9.py"
  --input-csv "${INPUT_CSV}"
  --output-dir "${OUTPUT_DIR}"
  --output-csv "${OUTPUT_CSV}"
  --sem-count "${SEM_COUNT}"
  --aging-hours "${AGING_HOURS}"
  --agglomeration-rate "${AGGLOMERATION_RATE}"
  --geo-sigma "${GEO_SIGMA}"
)

if [[ -n "${MODEL_JSON}" ]]; then
  CMD+=(--model "${MODEL_JSON}")
fi

"${CMD[@]}"

echo

echo "Done. Predicted table: ${OUTPUT_CSV}"
echo "Per-run JSON files: ${OUTPUT_DIR}"
