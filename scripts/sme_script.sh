#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

LAMMPS_BIN="${LAMMPS_BIN:-/usr/local/bin/lmp}"

MINIMIZE_IN="${MINIMIZE_IN:-${PROJECT_ROOT}/src/lmps/minimize.lmp}"
THERMAL_IN="${THERMAL_IN:-${PROJECT_ROOT}/src/lmps/thermal_cycle.lmp}"
SE_IN="${SE_IN:-${PROJECT_ROOT}/src/lmps/se.lmp}"
SME_IN="${SME_IN:-${PROJECT_ROOT}/src/lmps/sme.lmp}"

POT_LIB="${POT_LIB:-${PROJECT_ROOT}/potentials/library.meam}"
POT_FILE="${POT_FILE:-${PROJECT_ROOT}/potentials/NiTi.meam}"

SEED="${SEED:-42}"

DATA_FILES=(
  "${PROJECT_ROOT}/data/structures/niti_2nm_Ni50_amorphous.data"
  # add your 5/10/15/20 nm data files here
)

RAW_OUT_ROOT="${RAW_OUT_ROOT:-${PROJECT_ROOT}/data/raw_output}"
LOG_DIR="${LOG_DIR:-${PROJECT_ROOT}/logs}"

THERMAL_LOW_RESTART_NAME="${THERMAL_LOW_RESTART_NAME:-restart.stage2_end}"
THERMAL_HIGH_RESTART_NAME="${THERMAL_HIGH_RESTART_NAME:-restart.stage4_end}"

mkdir -p "${LOG_DIR}"
mkdir -p "${RAW_OUT_ROOT}"

[[ -x "${LAMMPS_BIN}" ]] || { echo "ERROR: LAMMPS not found: ${LAMMPS_BIN}" >&2; exit 1; }
[[ -f "${POT_LIB}" && -f "${POT_FILE}" ]] || { echo "ERROR: Missing potential files" >&2; exit 1; }
[[ -f "${MINIMIZE_IN}" ]] || { echo "ERROR: Missing ${MINIMIZE_IN}" >&2; exit 1; }
[[ -f "${THERMAL_IN}" ]]   || { echo "ERROR: Missing ${THERMAL_IN}" >&2; exit 1; }
[[ -f "${SE_IN}" ]]        || { echo "ERROR: Missing ${SE_IN}" >&2; exit 1; }
[[ -f "${SME_IN}" ]]       || { echo "ERROR: Missing ${SME_IN}" >&2; exit 1; }

get_base_name() { basename "$1" .data; }

run_lammps() {
  local logfile=$1; shift
  mpirun -n "${MPI_RANKS}" "${LAMMPS_BIN}" "$@" 2>&1 | tee "${logfile}"
}

TOTAL_CORES="$(nproc)"
if [[ ${TOTAL_CORES} -ge 16 ]]; then
  export OMP_NUM_THREADS=2; MPI_RANKS=8
elif [[ ${TOTAL_CORES} -ge 8 ]]; then
  export OMP_NUM_THREADS=2; MPI_RANKS=4
elif [[ ${TOTAL_CORES} -ge 4 ]]; then
  export OMP_NUM_THREADS=2; MPI_RANKS=2
else
  export OMP_NUM_THREADS=1; MPI_RANKS=${TOTAL_CORES}
fi

echo "INFO: cores=${TOTAL_CORES} MPI_RANKS=${MPI_RANKS} OMP_NUM_THREADS=${OMP_NUM_THREADS}"

# ----------------------------
# 1) Minimization
# ----------------------------
echo ""
echo "[1/4] Minimization (if missing)..."

for datafile in "${DATA_FILES[@]}"; do
  [[ -f "${datafile}" ]] || { echo "WARN: missing ${datafile}"; continue; }

  name="$(get_base_name "${datafile}")"
  min_dir="${RAW_OUT_ROOT}/${name}/minimize"
  min_restart_dir="${min_dir}/restarts"
  minimized="${min_dir}/minimized.data"
  marker="${min_dir}/.minimize_complete"

  if [[ -f "${minimized}" && -f "${marker}" ]]; then
    echo "  ✓ ${name} minimized exists"
    continue
  fi

  mkdir -p "${min_dir}" "${min_restart_dir}"
  logfile="${LOG_DIR}/minimize_${name}.log"

  echo "  → minimizing ${name}"
  run_lammps "${logfile}" \
    -var restart_file "none" \
    -var input_data "${datafile}" \
    -var output_dir "${min_dir}" \
    -var restart_dir "${min_restart_dir}" \
    -var pot_lib "${POT_LIB}" \
    -var pot_file "${POT_FILE}" \
    -var seed "${SEED}" \
    -in "${MINIMIZE_IN}"

  [[ -f "${minimized}" ]] || { echo "ERROR: minimize did not create ${minimized}" >&2; exit 1; }
  touch "${marker}"
done

# ----------------------------
# 2) Thermal cycle
# ----------------------------
echo ""
echo "[2/4] Thermal cycle..."

for datafile in "${DATA_FILES[@]}"; do
  [[ -f "${datafile}" ]] || continue

  name="$(get_base_name "${datafile}")"
  minimized="${RAW_OUT_ROOT}/${name}/minimize/minimized.data"
  thermal_dir="${RAW_OUT_ROOT}/${name}/thermal"
  thermal_restart_dir="${thermal_dir}/restarts"
  marker="${thermal_dir}/.thermal_complete"
  check_file="${thermal_dir}/${THERMAL_HIGH_RESTART_NAME}"

  if [[ -f "${marker}" ]]; then
    echo "  ✓ ${name} thermal already done"
    continue
  fi

  [[ -f "${minimized}" ]] || { echo "ERROR: missing ${minimized}" >&2; exit 1; }

  mkdir -p "${thermal_dir}" "${thermal_restart_dir}"
  logfile="${LOG_DIR}/thermal_${name}.log"

  echo "  → thermal ${name}"
  run_lammps "${logfile}" \
    -var restart_file "none" \
    -var input_data "${minimized}" \
    -var output_dir "${thermal_dir}" \
    -var restart_dir "${thermal_restart_dir}" \
    -var pot_lib "${POT_LIB}" \
    -var pot_file "${POT_FILE}" \
    -var seed "${SEED}" \
    -in "${THERMAL_IN}"

  [[ -f "${check_file}" ]] || { echo "ERROR: thermal did not create ${check_file}" >&2; exit 1; }
  touch "${marker}"
done

# ----------------------------
# 3) SE from high-T restart
# ----------------------------
echo ""
echo "[3/4] SE (from high-T thermal restart)..."

for datafile in "${DATA_FILES[@]}"; do
  [[ -f "${datafile}" ]] || continue

  name="$(get_base_name "${datafile}")"
  thermal_dir="${RAW_OUT_ROOT}/${name}/thermal"
  restart_hi="${thermal_dir}/${THERMAL_HIGH_RESTART_NAME}"

  se_dir="${RAW_OUT_ROOT}/${name}/se"
  se_restart_dir="${se_dir}/restarts"
  marker="${se_dir}/.se_complete"

  if [[ -f "${marker}" ]]; then
    echo "  ✓ ${name} SE already done"
    continue
  fi

  [[ -f "${restart_hi}" ]] || { echo "ERROR: missing ${restart_hi}" >&2; exit 1; }

  mkdir -p "${se_dir}" "${se_restart_dir}"
  logfile="${LOG_DIR}/se_${name}.log"

  echo "  → SE ${name}"
  run_lammps "${logfile}" \
    -var restart_file "${restart_hi}" \
    -var input_data "unused.data" \
    -var output_dir "${se_dir}" \
    -var restart_dir "${se_restart_dir}" \
    -var pot_lib "${POT_LIB}" \
    -var pot_file "${POT_FILE}" \
    -var seed "${SEED}" \
    -in "${SE_IN}"

  touch "${marker}"
done

# ----------------------------
# 4) SME from low-T restart
# ----------------------------
echo ""
echo "[4/4] SME (from low-T thermal restart)..."

for datafile in "${DATA_FILES[@]}"; do
  [[ -f "${datafile}" ]] || continue

  name="$(get_base_name "${datafile}")"
  thermal_dir="${RAW_OUT_ROOT}/${name}/thermal"
  restart_lo="${thermal_dir}/${THERMAL_LOW_RESTART_NAME}"

  sme_dir="${RAW_OUT_ROOT}/${name}/sme"
  sme_restart_dir="${sme_dir}/restarts"
  marker="${sme_dir}/.sme_complete"

  if [[ -f "${marker}" ]]; then
    echo "  ✓ ${name} SME already done"
    continue
  fi

  [[ -f "${restart_lo}" ]] || { echo "ERROR: missing ${restart_lo}" >&2; exit 1; }

  mkdir -p "${sme_dir}" "${sme_restart_dir}"
  logfile="${LOG_DIR}/sme_${name}.log"

  echo "  → SME ${name}"
  run_lammps "${logfile}" \
    -var restart_file "${restart_lo}" \
    -var input_data "unused.data" \
    -var output_dir "${sme_dir}" \
    -var restart_dir "${sme_restart_dir}" \
    -var pot_lib "${POT_LIB}" \
    -var pot_file "${POT_FILE}" \
    -var seed "${SEED}" \
    -in "${SME_IN}"

  touch "${marker}"
done

echo ""
echo "=========================================="
echo "✓ COMPLETE: minimize + thermal + SE + SME"
echo "Logs: ${LOG_DIR}"
echo "Output: ${RAW_OUT_ROOT}"
echo "=========================================="