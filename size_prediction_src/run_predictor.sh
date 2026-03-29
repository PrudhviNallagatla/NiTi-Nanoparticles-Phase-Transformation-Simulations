#!/usr/bin/env bash
set -euo pipefail

# Simple one-file runner for the NiTi size predictor.
# Edit values below, save, then run: ./run_predictor.sh

PYTHON_BIN="python3"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREDICTOR_PY="${SCRIPT_DIR}/predictor.py"

# ------------------------------
# Editable process parameters
# ------------------------------
VOLTAGE_V=25
CURRENT_A=3
DUTY_PERCENT=30
FREQ_KHZ=5

# Optional observation parameters
SEM_COUNT=10
AGING_HOURS=0.0
AGGLOMERATION_RATE=0.08
GEO_SIGMA=1.35

# Optional model I/O
MODEL_JSON=""
SAVE_MODEL_JSON=""
OUTPUT_JSON="${SCRIPT_DIR}/prediction_output.json"

CMD=(
  "${PYTHON_BIN}" "${PREDICTOR_PY}"
  --voltage "${VOLTAGE_V}"
  --current "${CURRENT_A}"
  --duty-percent "${DUTY_PERCENT}"
  --freq-khz "${FREQ_KHZ}"
  --sem-count "${SEM_COUNT}"
  --aging-hours "${AGING_HOURS}"
  --agglomeration-rate "${AGGLOMERATION_RATE}"
  --geo-sigma "${GEO_SIGMA}"
)

if [[ -n "${MODEL_JSON}" ]]; then
  CMD+=(--model "${MODEL_JSON}")
fi

if [[ -n "${SAVE_MODEL_JSON}" ]]; then
  CMD+=(--save-model "${SAVE_MODEL_JSON}")
fi

mkdir -p "$(dirname "${OUTPUT_JSON}")"
"${CMD[@]}" | tee "${OUTPUT_JSON}"

echo

echo "Prediction written to: ${OUTPUT_JSON}"
