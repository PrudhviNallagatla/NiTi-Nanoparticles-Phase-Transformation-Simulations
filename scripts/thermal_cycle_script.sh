#!/bin/bash
# ==========================================================================
# AUTOMATED THERMAL CYCLING WORKFLOW WITH CRASH RECOVERY
# Purpose: Complete pipeline from raw structure → minimization → thermal characterization
# Features: Auto-restart from checkpoints, dynamic resource allocation, progress tracking
# ==========================================================================

set -euo pipefail  # Exit on error, undefined variables, or pipe failures

# ============================================
# 1. PROJECT CONFIGURATION
# ============================================

# Input structures to process (add/remove as needed)
DATA_FILES=(
    "/home/rimuru/workspace/data/structures/niti_10nm_Ni50_single_crystal.data"
)

PROJECT_ROOT="/home/rimuru/workspace"
LAMMPS_BIN="/usr/local/bin/lmp"
SEED=42  # Random seed for reproducible velocity initialization

# MEAM Potential Files (Ni-Ti binary system)
POT_LIB="/home/rimuru/workspace/potentials/library.meam"
POT_FILE="/home/rimuru/workspace/potentials/NiTi.meam"

# Hardware Detection (used for automatic parallelization)
TOTAL_CORES=$(nproc)
TOTAL_MEM_GB=$(awk '/MemTotal/ {printf "%.0f\n", $2/1024/1024}' /proc/meminfo)

# ============================================
# 2. DYNAMIC PERFORMANCE OPTIMIZATION
# Purpose: Balance MPI ranks vs OpenMP threads based on available hardware
# Strategy: MEAM potentials benefit from hybrid parallelization
# ============================================
echo "INFO: Detected ${TOTAL_CORES} CPU cores, ${TOTAL_MEM_GB} GB RAM"

if [[ ${TOTAL_CORES} -ge 16 ]]; then
    # High-end system: Hybrid MPI+OpenMP for best MEAM performance
    export OMP_NUM_THREADS=2
    MPI_RANKS=8
    echo "INFO: High-performance profile: 8 MPI ranks × 2 OpenMP threads"
    
elif [[ ${TOTAL_CORES} -eq 8 ]]; then
    # Standard workstation: Pure MPI (simpler, more stable)
    export OMP_NUM_THREADS=1
    MPI_RANKS=8
    echo "INFO: Standard profile: 8 MPI ranks × 1 OpenMP thread"
    
elif [[ ${TOTAL_CORES} -eq 4 ]]; then
    # Low-end system: Conservative hybrid approach
    export OMP_NUM_THREADS=2
    MPI_RANKS=2
    echo "INFO: Low-resource profile: 2 MPI ranks × 2 OpenMP threads"
    
else
    # Fallback: Use all cores as MPI ranks
    export OMP_NUM_THREADS=1
    MPI_RANKS=${TOTAL_CORES}
    echo "INFO: Fallback pure-MPI profile: ${TOTAL_CORES} MPI ranks"
fi

# ============================================
# 3. PRE-FLIGHT VALIDATION
# Purpose: Catch configuration errors before wasting compute time
# ============================================
if [[ ! -x "${LAMMPS_BIN}" ]]; then
    echo "ERROR: LAMMPS binary not found or not executable: ${LAMMPS_BIN}" >&2
    exit 1
fi

if [[ ! -f "${POT_LIB}" ]]; then
    echo "ERROR: MEAM library file missing: ${POT_LIB}" >&2
    exit 1
fi

if [[ ! -f "${POT_FILE}" ]]; then
    echo "ERROR: MEAM parameter file missing: ${POT_FILE}" >&2
    exit 1
fi

cd "${PROJECT_ROOT}"
mkdir -p logs

# ============================================
# 4. UTILITY FUNCTIONS
# ============================================

# Extract base filename without extension
get_base_name() {
    basename "$1" .data
}

# Parse atom count from LAMMPS data file header
get_atom_count() {
    local datafile=$1
    [[ ! -f "${datafile}" ]] && echo "0" && return
    local count
    count=$(grep -m1 "atoms" "${datafile}" 2>/dev/null | awk '{print $1}')
    echo "${count:-0}"
}

# Estimate memory requirement (rough heuristic: ~500 bytes/atom)
estimate_memory() {
    local atoms=$1
    # Returns memory in MB
    python3 -c "print(max(1, int(${atoms} * 500 / 1024 / 1024) + 1))"
}

# ============================================
# STEP 1: STRUCTURE ANALYSIS
# Purpose: Preview computational requirements before starting
# ============================================
echo ""
echo "========================================"
echo "[1/3] ANALYZING INPUT STRUCTURES"
echo "========================================"

for datafile in "${DATA_FILES[@]}"; do
    if [[ ! -f "${datafile}" ]]; then
        echo "WARNING: File not found (will skip): ${datafile}" >&2
        continue
    fi
    
    structure_name=$(get_base_name "${datafile}")
    atoms=$(get_atom_count "${datafile}")
    mem=$(estimate_memory "${atoms}")
    
    echo "Structure: ${structure_name}"
    echo "  ├─ Atoms:     ${atoms}"
    echo "  ├─ Est. RAM:  ~${mem} MB"
    echo "  └─ MPI Ranks: ${MPI_RANKS}"
    echo ""
done

# ============================================
# STEP 2: ENERGY MINIMIZATION
# Purpose: Relax atomic positions to stable configuration
# Strategy: Skip if already completed (idempotent design)
# ============================================
echo ""
echo "========================================"
echo "[2/3] ENERGY MINIMIZATION PHASE"
echo "========================================"
echo "Method: Multi-stage (SD → Anneal → FIRE → CG)"
echo "Target: Fmax < 0.01 meV/Å"
echo "========================================"
echo ""

for datafile in "${DATA_FILES[@]}"; do
    [[ ! -f "${datafile}" ]] && continue
    
    structure_name=$(get_base_name "${datafile}")
    
    # Output directory for this structure
    outdir="data/raw_output/${structure_name}/minimize"
    marker="${outdir}/.minimize_complete"  # Completion flag
    output_check="${outdir}/final_minimized.data"  # Expected output file
    
    # Check if already completed
    if [[ -f "${marker}" ]]; then
        echo "✓ ${structure_name} (already minimized, skipping)"
        continue
    fi
    
    mkdir -p "${outdir}"
    echo "→ Minimizing ${structure_name}..."
    logfile="logs/min_${structure_name}.log"

    # Execute LAMMPS minimization
    STRUCTURE_FILE="${datafile}" \
    OUTPUT_DIR="${outdir}" \
    POT_LIB="${POT_LIB}" \
    POT_FILE="${POT_FILE}" \
    mpirun -n "${MPI_RANKS}" "${LAMMPS_BIN}" \
        -in src/lmps/minimize.lmp \
        2>&1 | tee "${logfile}"
    
    # Validate output
    if [[ -f "${output_check}" ]]; then
        touch "${marker}"  # Mark as complete
        echo "  ✓ Minimization complete"
        echo ""
    else
        echo "  ✗ FAILED: Output file not created" >&2
        echo "  Check log: ${logfile}" >&2
        exit 1
    fi
done

# ============================================
# STEP 3: THERMAL CYCLING CHARACTERIZATION
# Purpose: Measure transformation temperatures (Ms, Mf, As, Af)
# Method: Cool → Heat cycle with continuous phase monitoring
# Features: Auto-resume from latest checkpoint on crash/restart
# ============================================
echo ""
echo "========================================"
echo "[3/3] THERMAL CYCLING CHARACTERIZATION"
echo "========================================"
echo "Temperature range: 200 K → 500 K → 200 K"
echo "Heating/Cooling rate: 0.35 K/ps"
echo "Output: Hysteresis curve (T vs Austenite Fraction)"
echo "========================================"
echo ""

for datafile in "${DATA_FILES[@]}"; do
    [[ ! -f "${datafile}" ]] && continue
    
    structure_name=$(get_base_name "${datafile}")

    # Define paths
    indata="data/raw_output/${structure_name}/minimize/final_minimized.data"
    outdir="data/raw_output/${structure_name}/thermal_cycle"
    restart_dir="${outdir}/restart_files"
    
    marker="${outdir}/.thermal_cycle_complete"  # Completion flag
    output_check="${restart_dir}/restart.stage4_end"  # Final stage checkpoint
    logfile="logs/thermal_${structure_name}.log"

    # Check if already completed
    if [[ -f "${marker}" ]]; then
        echo "✓ ${structure_name} (thermal cycle already complete)"
        continue
    fi

    # Validate input exists
    if [[ ! -f "${indata}" ]]; then
        echo "ERROR: Minimized structure not found: ${indata}" >&2
        echo "  Run minimization step first!" >&2
        exit 1
    fi

    mkdir -p "${outdir}" "${restart_dir}"

    # =====================================
    # SMART RESTART LOGIC
    # Purpose: Resume from most recent checkpoint to recover from crashes
    # Strategy: Look for stage_end (definitive) first, then rolling backups
    # =====================================
    
    # Find most recent restart file (sorted by modification time)
    latest_restart=$(ls -1t \
        "${restart_dir}"/restart.stage*_end \
        "${restart_dir}"/restart.rolling.* \
        2>/dev/null | head -n 1 || true)

    if [[ -n "${latest_restart}" ]]; then
        # Resume from checkpoint
        echo "→ ${structure_name}: Resuming from checkpoint"
        echo "  Restart file: $(basename "${latest_restart}")"
        START_CMD="read_restart ${latest_restart}"
    else
        # Start fresh
        echo "→ ${structure_name}: Starting thermal cycle from minimized structure"
        START_CMD="read_data ${indata}"
    fi

    # Execute LAMMPS thermal cycle
    OUTPUT_DIR="${outdir}" \
    RESTART_DIR="${restart_dir}" \
    START_CMD="${START_CMD}" \
    SEED="${SEED}" \
    POT_LIB="${POT_LIB}" \
    POT_FILE="${POT_FILE}" \
    mpirun -n "${MPI_RANKS}" "${LAMMPS_BIN}" \
        -in src/lmps/thermal_cycle.lmp \
        2>&1 | tee "${logfile}"

    # Validate completion
    if [[ -f "${output_check}" ]]; then
        touch "${marker}"  # Mark as complete
        echo "  ✓ Thermal cycle complete"
        echo "  Output: ${outdir}/hysteresis_full.txt"
        echo ""
    else
        echo "  ✗ FAILED: Simulation did not complete all stages" >&2
        echo "  Check log: ${logfile}" >&2
        echo "  To retry: Re-run this script (will auto-resume)" >&2
        exit 1
    fi
done

echo ""
echo "========================================"
echo "✓ ALL SIMULATIONS COMPLETE"
echo "========================================"
echo "Next steps:"
echo "1. Analyze hysteresis curves: data/raw_output/*/thermal_cycle/hysteresis_full.txt"
echo "2. Extract transformation temperatures (Ms, Mf, As, Af)"
echo "3. Visualize trajectories: ovito data/raw_output/*/thermal_cycle/dump.cycle.*.gz"
echo "========================================"