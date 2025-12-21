#!/usr/bin/env bash
set -euo pipefail

# --- CONFIG ---
PROJECT_ROOT="/home/rimuru/workspace"
LAMMPS_BIN="/usr/local/bin/lmp"

# Inputs
MINIMIZE_IN="${PROJECT_ROOT}/src/lmps/minimize.lmp"
THERMAL_IN="${PROJECT_ROOT}/src/lmps/thermal_cycle.lmp"
SE_IN="${PROJECT_ROOT}/src/lmps/se.lmp"
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
# FORCE 8 CORES (As requested)
export MPI_RANKS=8

echo ">>> RUNNING WITH ${MPI_RANKS} CORES"

# ----------------------------
# 1) MINIMIZATION
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
# 2) THERMAL
# ----------------------------
for datafile in "${DATA_FILES[@]}"; do
  name=$(basename "${datafile}" .data)
  dir="${RAW_OUT_ROOT}/${name}/thermal_cycle"
  mkdir -p "${dir}/restart_files"
  
  if [[ -f "${dir}/restart_files/restart.stage4_end" ]]; then
    echo ">>> Skipping Thermal (Found)"
  else
    echo ">>> Thermal Cycle ${name}"
    INPUT_DATA="${RAW_OUT_ROOT}/${name}/minimize/final_minimized.data" \
    OUTPUT_DIR="${dir}" RESTART_DIR="${dir}/restart_files" \
    POT_LIB="${POT_LIB}" POT_FILE="${POT_FILE}" SEED="42" \
    mpirun -n "${MPI_RANKS}" "${LAMMPS_BIN}" -in "${THERMAL_IN}" 2>&1 | tee "${LOG_DIR}/thermal_${name}.log"
  fi
done

# ----------------------------
# 3) SME (Shape Memory Effect) - PRIORITY
# ----------------------------
for datafile in "${DATA_FILES[@]}"; do
  name=$(basename "${datafile}" .data)
  dir="${RAW_OUT_ROOT}/${name}/sme"
  # SME starts from Low Temp (Stage 2 End)
  restart_lo="${RAW_OUT_ROOT}/${name}/thermal_cycle/restart_files/restart.stage2_end"
  
  mkdir -p "${dir}/restarts"
  
  if [[ -f "${dir}/final.data" ]]; then
    echo ">>> Skipping SME (Found)"
  else
    echo ">>> Running SME ${name}"
    # VISIBLE OUTPUT: Uses | tee
    INPUT_DATA="unused" OUTPUT_DIR="${dir}" RESTART_DIR="${dir}/restarts" \
    START_CMD="read_restart ${restart_lo}" \
    POT_LIB="${POT_LIB}" POT_FILE="${POT_FILE}" SEED="42" \
    mpirun -n "${MPI_RANKS}" "${LAMMPS_BIN}" -in "${SME_IN}" 2>&1 | tee "${LOG_DIR}/sme_${name}.log"
  fi
done

# ----------------------------
# 4) SE (Superelasticity) - AFTER SME
# ----------------------------
for datafile in "${DATA_FILES[@]}"; do
  name=$(basename "${datafile}" .data)
  dir="${RAW_OUT_ROOT}/${name}/se"
  # SE starts from High Temp (Stage 4 End)
  restart_hi="${RAW_OUT_ROOT}/${name}/thermal_cycle/restart_files/restart.stage4_end"
  
  mkdir -p "${dir}/restarts"
  
  if [[ -f "${dir}/final.data" ]]; then
    echo ">>> Skipping SE (Found)"
  else
    echo ">>> Running SE ${name}"
    # VISIBLE OUTPUT: Uses | tee to show you the progress
    INPUT_DATA="unused" OUTPUT_DIR="${dir}" RESTART_DIR="${dir}/restarts" \
    THERMAL_RESTART="${restart_hi}" \
    POT_LIB="${POT_LIB}" POT_FILE="${POT_FILE}" SEED="42" \
    mpirun -n "${MPI_RANKS}" "${LAMMPS_BIN}" -in "${SE_IN}" 2>&1 | tee "${LOG_DIR}/se_${name}.log"
  fi
done

echo ">>> ALL DONE."