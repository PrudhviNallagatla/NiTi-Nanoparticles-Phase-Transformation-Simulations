#!/usr/bin/env bash
set -euo pipefail

# Linux runner for size_predictor/predictor.py
# Example:
#   ./scripts/run_size_predictor.sh --voltage 120 --capacitance 1e-7 --duty 0.35 --pulse-on-us 8 --freq 40000

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"

usage() {
  cat <<'EOF'
Usage:
  run_size_predictor.sh --voltage V [legacy or experimental args] [options]

Required:
  --voltage <float>         Open-gap voltage (V)

Legacy mode args:
  --capacitance <float>     Capacitance (F)
  --duty <float>            Duty cycle in (0,1)
  --pulse-on-us <float>     Pulse-on time (microseconds)
  --freq <float>            Pulse frequency (Hz), default 50000

Experimental mode args:
  --current <float>         Current (A)
  --duty-percent <float>    Duty cycle (%)
  --freq-khz <float>        Frequency (kHz)
  --capacitance <float>     Optional: if omitted, equivalent capacitance is inferred

Optional:
  --sem-count <int>         Number of particles identified in SEM, default 10
  --aging-hours <float>     Delay between synthesis and SEM in hours, default 0
  --agglomeration-rate <f>  Effective agglomeration growth per hour, default 0.08
  --geo-sigma <float>       Geometric sigma for size spread, default 1.35
  --model <path>            Existing model JSON to load
  --calibrate-csv <path>    CSV with measured data for calibration
  --save-model <path>       Save resulting model JSON
  --output <path>           Save output JSON (default: stdout only)
  --python <path>           Python executable (overrides PYTHON_BIN)
  -h, --help                Show this help

Calibration CSV columns:
  legacy: voltage_v,capacitance_f,duty_cycle,pulse_on_us,measured_size_nm[,pulse_frequency_hz]
  experimental: voltage_v,current_a,duty_percent,frequency_khz,measured_size_nm[,capacitance_f]
EOF
}

VOLTAGE=""
CAPACITANCE=""
DUTY=""
PULSE_ON_US=""
FREQ="50000"
CURRENT=""
DUTY_PERCENT=""
FREQ_KHZ=""
SEM_COUNT="10"
AGING_HOURS="0"
AGGLOMERATION_RATE="0.08"
GEO_SIGMA="1.35"
MODEL=""
CAL_CSV=""
SAVE_MODEL=""
OUTPUT_FILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --voltage)
      VOLTAGE="$2"; shift 2 ;;
    --capacitance)
      CAPACITANCE="$2"; shift 2 ;;
    --duty)
      DUTY="$2"; shift 2 ;;
    --pulse-on-us)
      PULSE_ON_US="$2"; shift 2 ;;
    --freq)
      FREQ="$2"; shift 2 ;;
    --current)
      CURRENT="$2"; shift 2 ;;
    --duty-percent)
      DUTY_PERCENT="$2"; shift 2 ;;
    --freq-khz)
      FREQ_KHZ="$2"; shift 2 ;;
    --sem-count)
      SEM_COUNT="$2"; shift 2 ;;
    --aging-hours)
      AGING_HOURS="$2"; shift 2 ;;
    --agglomeration-rate)
      AGGLOMERATION_RATE="$2"; shift 2 ;;
    --geo-sigma)
      GEO_SIGMA="$2"; shift 2 ;;
    --model)
      MODEL="$2"; shift 2 ;;
    --calibrate-csv)
      CAL_CSV="$2"; shift 2 ;;
    --save-model)
      SAVE_MODEL="$2"; shift 2 ;;
    --output)
      OUTPUT_FILE="$2"; shift 2 ;;
    --python)
      PYTHON_BIN="$2"; shift 2 ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 2 ;;
  esac
done

if [[ -z "${VOLTAGE}" ]]; then
  echo "Missing required --voltage." >&2
  usage
  exit 2
fi

LEGACY_MODE=false
EXPERIMENTAL_MODE=false

if [[ -n "${CAPACITANCE}" && -n "${DUTY}" && -n "${PULSE_ON_US}" ]]; then
  LEGACY_MODE=true
fi

if [[ -n "${CURRENT}" && -n "${DUTY_PERCENT}" && -n "${FREQ_KHZ}" ]]; then
  EXPERIMENTAL_MODE=true
fi

if [[ "${LEGACY_MODE}" == false && "${EXPERIMENTAL_MODE}" == false ]]; then
  echo "Missing required arguments." >&2
  echo "Provide either legacy mode or experimental mode arguments." >&2
  usage
  exit 2
fi

PREDICTOR_PY="${PROJECT_ROOT}/size_predictor/predictor.py"
if [[ ! -f "${PREDICTOR_PY}" ]]; then
  echo "Predictor file not found at ${PREDICTOR_PY}" >&2
  exit 1
fi

CMD=(
  "${PYTHON_BIN}" "${PREDICTOR_PY}"
  --voltage "${VOLTAGE}"
  --sem-count "${SEM_COUNT}"
  --aging-hours "${AGING_HOURS}"
  --agglomeration-rate "${AGGLOMERATION_RATE}"
  --geo-sigma "${GEO_SIGMA}"
)

if [[ "${EXPERIMENTAL_MODE}" == true ]]; then
  CMD+=(
    --current "${CURRENT}"
    --duty-percent "${DUTY_PERCENT}"
    --freq-khz "${FREQ_KHZ}"
  )
  if [[ -n "${CAPACITANCE}" ]]; then
    CMD+=(--capacitance "${CAPACITANCE}")
  fi
else
  CMD+=(
    --capacitance "${CAPACITANCE}"
    --duty "${DUTY}"
    --pulse-on-us "${PULSE_ON_US}"
    --freq "${FREQ}"
  )
fi

if [[ -n "${MODEL}" ]]; then
  CMD+=(--model "${MODEL}")
fi

if [[ -n "${CAL_CSV}" ]]; then
  CMD+=(--calibrate-csv "${CAL_CSV}")
fi

if [[ -n "${SAVE_MODEL}" ]]; then
  CMD+=(--save-model "${SAVE_MODEL}")
fi

if [[ -n "${OUTPUT_FILE}" ]]; then
  mkdir -p "$(dirname "${OUTPUT_FILE}")"
  "${CMD[@]}" | tee "${OUTPUT_FILE}"
else
  "${CMD[@]}"
fi
