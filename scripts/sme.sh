#!/usr/bin/env bash
set -euo pipefail

# --- CONFIG ---
PROJECT_ROOT="/home/rimuru/workspace"
LAMMPS_BIN="/usr/local/bin/lmp"

# Inputs
MINIMIZE_IN="${PROJECT_ROOT}/src/lmps/minimize.lmp"
THERMAL_IN="${PROJECT_ROOT}/src/lmps/thermal_cycle.lmp"
SME_IN="${PROJECT_ROOT}/src/lmps/sme.lmp"

POT_LIB="${PROJECT_ROOT}/potentials/library.meam"
POT_FILE="${PROJECT_ROOT}/potentials/NiTi.meam"

DATA_FILES=("${PROJECT_ROOT}/data/structures/niti_10nm_Ni50_single_crystal.data")

# Outputs
RAW_OUT_ROOT="${PROJECT_ROOT}/data/raw_output"
LOG_DIR="${PROJECT_ROOT}/logs"
mkdir -p "${LOG_DIR}" "${RAW_OUT_ROOT}"

# Resources
export OMP_NUM_THREADS=1
export MPI_RANKS=8

echo ">>> RUNNING SME WORKFLOW WITH ${MPI_RANKS} CORES"

# ----------------------------
# 1) MINIMIZATION (Prerequisite for Thermal)
# ----------------------------
for datafile in "${DATA_FILES[@]}"; do
  name=$(basename "${datafile}" .data)
  dir="${RAW_OUT_ROOT}/${name}/minimize"
  mkdir -p "${dir}/restarts"
  
  if [[ -f "${dir}/final_minimized.data" ]]; then
    echo ">>> Skipping Min (Found)"
  else
    echo ">>> Minimizing ${name}"
    STRUCTURE_FILE="${datafile}" OUTPUT_DIR="${dir}" POT_LIB="${POT_LIB}" POT_FILE="${POT_FILE}" \
    mpirun -n "${MPI_RANKS}" "${LAMMPS_BIN}" -in "${MINIMIZE_IN}" 2>&1 | tee "${LOG_DIR}/minimize_${name}.log"
  fi
done

# ----------------------------
# 2) THERMAL (Prerequisite for SME - Generates Low T State)
# ----------------------------
for datafile in "${DATA_FILES[@]}"; do
  name=$(basename "${datafile}" .data)
  dir="${RAW_OUT_ROOT}/${name}/thermal_cycle"
  mkdir -p "${dir}/restart_files"
  
  if [[ -f "${dir}/restart_files/restart.stage2_end" ]]; then
    echo ">>> Skipping Thermal (Low Temp State Found)"
  else
    echo ">>> Thermal Cycle ${name}"
    INPUT_DATA="${RAW_OUT_ROOT}/${name}/minimize/final_minimized.data" \
    OUTPUT_DIR="${dir}" RESTART_DIR="${dir}/restart_files" \
    POT_LIB="${POT_LIB}" POT_FILE="${POT_FILE}" SEED="42" \
    mpirun -n "${MPI_RANKS}" "${LAMMPS_BIN}" -in "${THERMAL_IN}" 2>&1 | tee "${LOG_DIR}/thermal_${name}.log"
  fi
done

# ----------------------------
# 3) SME (Shape Memory Effect)
# ----------------------------

for datafile in "${DATA_FILES[@]}"; do
  name=$(basename "${datafile}" .data)
  dir="${RAW_OUT_ROOT}/${name}/sme"
  # SME starts from Low Temp (Stage 2 End of Thermal Cycle)
  restart_lo="${RAW_OUT_ROOT}/${name}/thermal_cycle/restart_files/restart.stage2_end"
  
  mkdir -p "${dir}/restarts"
  
  if [[ ! -f "${restart_lo}" ]]; then
    echo "CRITICAL ERROR: Thermal restart file not found for ${name}. Cannot run SME."
    exit 1
  fi

  if [[ -f "${dir}/final.data" ]]; then
    echo ">>> Skipping SME (Found)"
  else
    echo ">>> Running SME ${name}"
    INPUT_DATA="unused" OUTPUT_DIR="${dir}" RESTART_DIR="${dir}/restarts" \
    START_CMD="read_restart ${restart_lo}" \
    POT_LIB="${POT_LIB}" POT_FILE="${POT_FILE}" SEED="42" \
    mpirun -n "${MPI_RANKS}" "${LAMMPS_BIN}" -in "${SME_IN}" 2>&1 | tee "${LOG_DIR}/sme_${name}.log"
  fi
done

echo ">>> SME WORKFLOW COMPLETED."