#!/bin/bash
# ==========================================================================
# ROBUST CHARACTERIZATION RUNNER (IQ 180 VERSION)
# ==========================================================================
# PURPOSE: Automates minimization and thermal cycling with crash recovery.
# DEPENDENCIES: 
#   - src/lmps/minimize.lmp (Standard Minimization)
#   - src/lmps/thermal_cycle.lmp (The Robust Cycle we just created)
# ==========================================================================

set -euo pipefail

# ============================================
# 1. CONFIGURATION
# ============================================
DATA_FILES=(
    "/home/rimuru/workspace/data/structures/niti_10nm_Ni50_single_crystal.data"
)

PROJECT_ROOT="/home/rimuru/workspace"
LAMMPS_BIN="/usr/local/bin/lmp"
SEED=42

# Potential Files
POT_LIB="/home/rimuru/workspace/potentials/library.meam"
POT_FILE="/home/rimuru/workspace/potentials/NiTi.meam"

# Hardware Detection
TOTAL_CORES=$(nproc)
TOTAL_MEM_GB=$(awk '/MemTotal/ {printf "%.0f\n", $2/1024/1024}' /proc/meminfo)

# ============================================
# 2. DYNAMIC PERFORMANCE SETTINGS
# ============================================
# Logic: Balance MPI ranks vs OpenMP threads based on core count
echo "INFO: Total cores detected: ${TOTAL_CORES}"

if [[ ${TOTAL_CORES} -ge 16 ]]; then
    # 16+ Cores: 8 MPI Ranks x 2 OpenMP Threads (Good for MEAM)
    export OMP_NUM_THREADS=2
    MPI_RANKS=8
    echo "INFO: Using 16+ core profile (MPI_RANKS=8, OMP_NUM_THREADS=2)"
    
elif [[ ${TOTAL_CORES} -eq 8 ]]; then
    # 8 Cores: 8 MPI Ranks x 1 OpenMP Thread (Standard)
    export OMP_NUM_THREADS=1
    MPI_RANKS=8
    echo "INFO: Using 8-core profile (MPI_RANKS=8, OMP_NUM_THREADS=1)"
    
elif [[ ${TOTAL_CORES} -eq 4 ]]; then
    # 4 Cores: 2 MPI Ranks x 2 OpenMP Threads (Safe fallback)
    export OMP_NUM_THREADS=2
    MPI_RANKS=2
    echo "INFO: Using 4-core profile (MPI_RANKS=2, OMP_NUM_THREADS=2)"
    
else
    # Fallback: Pure MPI
    export OMP_NUM_THREADS=1
    MPI_RANKS=${TOTAL_CORES}
    echo "INFO: Using fallback PURE MPI profile (MPI_RANKS=${TOTAL_CORES})"
fi

# ============================================
# 3. PRE-FLIGHT CHECKS
# ============================================
if [[ ! -x "${LAMMPS_BIN}" ]]; then
    echo "ERROR: LAMMPS binary not found at ${LAMMPS_BIN}" >&2; exit 1
fi
if [[ ! -f "${POT_LIB}" || ! -f "${POT_FILE}" ]]; then
    echo "ERROR: Potential files not found." >&2; exit 1
fi

cd "${PROJECT_ROOT}"
mkdir -p logs

# ============================================
# 4. HELPER FUNCTIONS
# ============================================
get_base_name() {
    basename "$1" .data
}

get_atom_count() {
    local datafile=$1
    [[ ! -f "${datafile}" ]] && echo "0" && return
    local count
    count=$(grep -m1 "atoms" "${datafile}" 2>/dev/null | awk '{print $1}')
    echo "${count:-0}"
}

estimate_memory() {
    local atoms=$1
    # Returns MB (integer) - simple safety buffer estimation
    python3 -c "print(max(1, int(${atoms} * 500 / 1024 / 1024) + 1))"
}

# ============================================
# STEP 1: ANALYZE STRUCTURES
# ============================================
echo ""
echo "[1/3] Analyzing structures..."
for datafile in "${DATA_FILES[@]}"; do
    if [[ ! -f "${datafile}" ]]; then
        echo "WARNING: File not found: ${datafile}" >&2; continue
    fi
    structure_name=$(get_base_name "${datafile}")
    atoms=$(get_atom_count "${datafile}")
    mem=$(estimate_memory "${atoms}")
    echo "   • ${structure_name}"
    echo "     Atoms: ${atoms} | RAM: ~${mem} MB | MPI Ranks: ${MPI_RANKS}"
done

# ============================================
# STEP 2: MINIMIZATION
# ============================================
echo ""
echo "[2/3] Minimizing structures (sequentially)..."
for datafile in "${DATA_FILES[@]}"; do
    [[ ! -f "${datafile}" ]] && continue
    structure_name=$(get_base_name "${datafile}")
    
    outdir="data/raw_output/${structure_name}/minimize"
    marker="${outdir}/.minimize_complete"
    output_check="${outdir}/minimized.data"
    
    if [[ -f "${marker}" ]]; then
        echo "   ✓ ${structure_name} (already done)"
        continue
    fi
    
    mkdir -p "${outdir}"
    echo "   → Running ${structure_name}..."
    logfile="logs/min_${structure_name}.log"

    # Run LAMMPS Minimize
    # NOTE: Ensure src/lmps/minimize.lmp exists and produces 'minimized.data'
    STRUCTURE_FILE="${datafile}" \
    OUTPUT_DIR="${outdir}" \
    POT_LIB="${POT_LIB}" \
    POT_FILE="${POT_FILE}" \
    mpirun -n "${MPI_RANKS}" "${LAMMPS_BIN}" \
        -in src/lmps/minimize.lmp \
        2>&1 | tee "${logfile}"
    
    if [[ -f "${output_check}" ]]; then
        touch "${marker}"
        echo "           ✓ ${structure_name} complete."
    else
        echo "           ✗ ${structure_name} FAILED (check ${logfile})" >&2
        exit 1
    fi
done

# ============================================
# STEP 3: THERMAL CYCLING (The Characterization)
# ============================================
echo ""
echo "[3/3] Running Thermal Cycle (sequentially)..."
for datafile in "${DATA_FILES[@]}"; do
    [[ ! -f "${datafile}" ]] && continue
    structure_name=$(get_base_name "${datafile}")

    # Define Input/Output paths
    indata="data/raw_output/${structure_name}/minimize/minimized.data"
    outdir="data/raw_output/${structure_name}/thermal_cycle"
    restart_dir="${outdir}/restart_files"
    
    marker="${outdir}/.thermal_cycle_complete"
    # We check for the final stage restart as proof of completion
    output_check="${restart_dir}/restart.stage4_end"
    logfile="logs/thermal_${structure_name}.log"

    if [[ -f "${marker}" ]]; then
        echo "   ✓ ${structure_name} (thermal cycle already done)"
        continue
    fi

    if [[ ! -f "${indata}" ]]; then
        echo "ERROR: Missing minimized input data: ${indata}" >&2
        exit 1
    fi

    mkdir -p "${outdir}" "${restart_dir}"

    # --- ROBUST RESTART LOGIC ---
    # Look for the newest restart file (stage_end or rolling)
    latest_restart=$(ls -1t "${restart_dir}"/restart.stage*_end "${restart_dir}"/restart.rolling.* 2>/dev/null | head -n 1 || true)

    if [[ -n "${latest_restart}" ]]; then
        echo "   ↳ Found restart checkpoint: $(basename "${latest_restart}")"
        START_CMD="read_restart ${latest_restart}"
    else
        echo "   ↳ Starting fresh from minimized structure..."
        START_CMD="read_data ${indata}"
    fi

    # Run LAMMPS Thermal Cycle
    OUTPUT_DIR="${outdir}" \
    RESTART_DIR="${restart_dir}" \
    START_CMD="${START_CMD}" \
    SEED="${SEED}" \
    POT_LIB="${POT_LIB}" \
    POT_FILE="${POT_FILE}" \
    mpirun -n "${MPI_RANKS}" "${LAMMPS_BIN}" \
        -in src/lmps/thermal_cycle.lmp \
        2>&1 | tee "${logfile}"

    if [[ -f "${output_check}" ]]; then
        touch "${marker}"
        echo "           ✓ ${structure_name} thermal cycle complete."
    else
        echo "           ✗ ${structure_name} FAILED (check ${logfile})" >&2
        exit 1
    fi
done

# ============================================
# SUMMARY & NEXT STEPS
# ============================================
echo ""
echo "=========================================="
echo "✓ PART 1 (CHARACTERIZATION) COMPLETE"
echo "=========================================="
echo "ACTION REQUIRED:"
echo "1. Plot the Hysteresis Loop from:"
echo "   data/raw_output/<structure_name>/thermal_cycle/hysteresis_full.txt"
echo "   (Columns: 4=Temp_global, 5=Temp_core, 6=PE, 13=Frac_Aus)"
echo ""
echo "2. Identify transformation temperatures (Mf, Af)."
echo "3. Update 'src/lmps/sme_deformation.lmp' with correct temperatures."
echo "4. Run ./scripts/run_02_deform.sh"
echo "=========================================="