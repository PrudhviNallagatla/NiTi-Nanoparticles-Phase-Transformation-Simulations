#!/bin/bash
# filepath: /home/rimuru/workspace/scripts/sme_scripts/run_01_characterize.sh
# ==========================================================================
# ROBUST CHARACTERIZATION RUNNER (IQ 180 VERSION)
# ==========================================================================
# MODIFIED: 
#  1. Matches new LAMMPS restart naming (stage*_end + rolling).
#  2. Safer memory estimation for MEAM/Tiled comms.
#  3. Auto-detects core count for hybrid MPI/OMP.
# ==========================================================================

set -euo pipefail

# ============================================
# CONFIGURATION
# ============================================
DATA_FILES=(
    "/home/rimuru/workspace/data/structures/niti_10nm_Ni50_single_crystal.data"
)

PROJECT_ROOT="/home/rimuru/workspace"
LAMMPS_BIN="/usr/local/bin/lmp"
SEED=42

POT_LIB="/home/rimuru/workspace/potentials/library.meam"
POT_FILE="/home/rimuru/workspace/potentials/NiTi.meam"

TOTAL_CORES=$(nproc)
TOTAL_MEM_GB=$(awk '/MemTotal/ {printf "%.0f\n", $2/1024/1024}' /proc/meminfo)

# ============================================
# 🚀 DYNAMIC PERFORMANCE SETTINGS
# ============================================
echo "INFO: Total cores detected: ${TOTAL_CORES}"

if [[ ${TOTAL_CORES} -ge 16 ]]; then
    # 16-core (or larger) VM: Use 8x2 hybrid
    export OMP_NUM_THREADS=2
    MPI_RANKS=8
    echo "INFO: Using 16+ core profile (MPI_RANKS=8, OMP_NUM_THREADS=2)"
    
elif [[ ${TOTAL_CORES} -eq 8 ]]; then
    # 8-core VM: Use 4x2 hybrid
    export OMP_NUM_THREADS=1
    MPI_RANKS=8
    echo "INFO: Using 8-core profile (MPI_RANKS=8, OMP_NUM_THREADS=1)"
    
elif [[ ${TOTAL_CORES} -eq 4 ]]; then
    # 4-core VM: Use 2x2 hybrid
    export OMP_NUM_THREADS=2
    MPI_RANKS=2
    echo "INFO: Using 4-core profile (MPI_RANKS=2, OMP_NUM_THREADS=2)"
    
else
    # Fallback for small/unusual VMs: Pure MPI (safe mode)
    export OMP_NUM_THREADS=1
    MPI_RANKS=${TOTAL_CORES}
    echo "INFO: Using fallback PURE MPI profile (MPI_RANKS=${TOTAL_CORES}, OMP_NUM_THREADS=1)"
fi

# ============================================
# Pre-flight Checks
# ============================================
if [[ ! -x "${LAMMPS_BIN}" ]]; then
    echo "ERROR: LAMMPS binary not found at ${LAMMPS_BIN}" >&2; exit 1
fi
if [[ ! -f "${POT_LIB}" || ! -f "${POT_FILE}" ]]; then
    echo "ERROR: Potential files not found." >&2; exit 1
fi

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "=========================================="
echo "SME Pipeline [PART 1]: CHARACTERIZATION"
echo "=========================================="
echo "Hardware detected:"
echo "   • CPU cores: ${TOTAL_CORES}"
echo "   • Total RAM: ${TOTAL_MEM_GB} GB"
echo "   • Files to process: ${#DATA_FILES[@]}"
echo "Performance:"
echo "   • OMP_NUM_THREADS: ${OMP_NUM_THREADS}"
echo "   • MPI_RANKS: ${MPI_RANKS}"
echo "=========================================="

# ============================================
# Helper Functions
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
    # MODIFIED: Increased to 500 bytes/atom for MEAM + Neighbor lists safety
    python3 -c "print(max(1, int(${atoms} * 500 / 1024 / 1024) + 1))"
}

# ============================================
# STEP 1: Analyze structures
# ============================================
echo ""
echo "[1/3] Analyzing structures..."
declare -A ATOM_COUNTS
declare -A MEM_ESTIMATES

for datafile in "${DATA_FILES[@]}"; do
    if [[ ! -f "${datafile}" ]]; then
        echo "WARNING: File not found: ${datafile}" >&2
        continue
    fi
    structure_name=$(get_base_name "${datafile}")
    atoms=$(get_atom_count "${datafile}")
    if [[ ${atoms} -eq 0 ]]; then
        echo "WARNING: Could not read atom count from ${datafile}" >&2
        continue
    fi
    mem=$(estimate_memory "${atoms}")
    ATOM_COUNTS["${structure_name}"]=${atoms}
    MEM_ESTIMATES["${structure_name}"]=${mem}
    echo "   • ${structure_name}"
    echo "       Atoms: ${atoms} | RAM: ~${mem} GB | MPI ranks: ${MPI_RANKS}"
done

# ============================================
# STEP 2: Minimize structures
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
    
    mpirun -n "${MPI_RANKS}" "${LAMMPS_BIN}" \
        -var structure_file "${datafile}" \
        -var output_dir "${outdir}" \
        -var pot_lib "${POT_LIB}" \
        -var pot_file "${POT_FILE}" \
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
echo "   ✓ All minimizations complete"

# ============================================
# STEP 3: Thermal Cycling (CHARACTERIZATION)
# ============================================
echo ""
echo "[3/3] Running Thermal Cycle (sequentially)..."
for datafile in "${DATA_FILES[@]}"; do
    [[ ! -f "${datafile}" ]] && continue
    structure_name=$(get_base_name "${datafile}")

    indata="data/raw_output/${structure_name}/minimize/minimized.data"
    outdir="data/raw_output/${structure_name}/thermal_cycle"
    restart_dir="${outdir}/restart_files"

    marker="${outdir}/.thermal_cycle_complete"
    output_check="${outdir}/final.data"
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
    # 1. Prefer completed stage boundaries (restart.stage*_end)
    # 2. Fallback to rolling checkpoints (restart.rolling.*)
    # 3. Default to none (start fresh)
    latest_restart=$(ls -1t "${restart_dir}"/restart.stage*_end "${restart_dir}"/restart.rolling.* 2>/dev/null | head -n 1 || true)
    
    if [[ -n "${latest_restart}" ]]; then
        echo "   ↳ Found restart checkpoint: $(basename "${latest_restart}")"
        RESTART_ARGS=(-var restart_file "${latest_restart}")
    else
        RESTART_ARGS=(-var restart_file none)
    fi

    mpirun -n "${MPI_RANKS}" "${LAMMPS_BIN}" \
        -var input_data "${indata}" \
        -var output_dir "${outdir}" \
        -var restart_dir "${restart_dir}" \
        -var seed ${SEED} \
        -var pot_lib "${POT_LIB}" \
        -var pot_file "${POT_FILE}" \
        "${RESTART_ARGS[@]}" \
        -in src/lmps/sme_thermal_cycle.lmp \
        2>&1 | tee "${logfile}"

    if [[ -f "${output_check}" ]]; then
        touch "${marker}"
        echo "           ✓ ${structure_name} thermal cycle complete."
    else
        echo "           ✗ ${structure_name} FAILED (check ${logfile})" >&2
        exit 1
    fi
done
echo "   ✓ Thermal cycling complete"

# ============================================
# Summary & NEXT STEPS
# ============================================
echo ""
echo "=========================================="
echo "✓ PART 1 (CHARACTERIZATION) COMPLETE"
echo "=========================================="
echo ""
echo "!!! ACTION REQUIRED !!!"
echo ""
echo "1.  Plot the Hysteresis Loop from:"
echo "    data/raw_output/<structure_name>/thermal_cycle/hysteresis.txt"
echo "    (Columns: 4=Temp_global, 5=Temp_core, 6=PE, 13=Frac_Aus)"
echo ""
echo "2.  Identify transformation temperatures (Mf, Af)."
echo ""
echo "3.  Update 'src/lmps/sme_deformation.lmp' with correct T_low/T_high."
echo ""
echo "4.  Run ./scripts/run_02_deform.sh"
echo "=========================================="