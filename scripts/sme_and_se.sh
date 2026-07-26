#!/usr/bin/env bash
set -euo pipefail

export OMP_NUM_THREADS=1
export MPI_RANKS=8

PRJ="/home/rimuru/workspace"
LMP="/usr/local/bin/lmp"
LIB="$PRJ/potentials/library.meam"
POT="$PRJ/potentials/NiTi.meam"

for df in "$PRJ/data/structures/niti_10nm_Ni50_single_crystal.data"; do
  name=$(basename "$df" .data)
  out="$PRJ/data/raw_output/$name"
  
  mkdir -p "$out/minimize/restarts" "$PRJ/logs"
  [[ -f "$out/minimize/final_minimized.data" ]] || STRUCTURE_FILE="$df" OUTPUT_DIR="$out/minimize" POT_LIB="$LIB" POT_FILE="$POT" mpirun -n $MPI_RANKS "$LMP" -in "$PRJ/src/lmps/minimize.lmp" 2>&1 | tee "$PRJ/logs/minimize_$name.log"

  mkdir -p "$out/thermal_cycle/restart_files"
  [[ -f "$out/thermal_cycle/restart_files/restart.stage4_end" ]] || INPUT_DATA="$out/minimize/final_minimized.data" OUTPUT_DIR="$out/thermal_cycle" RESTART_DIR="$out/thermal_cycle/restart_files" POT_LIB="$LIB" POT_FILE="$POT" SEED=42 mpirun -n $MPI_RANKS "$LMP" -in "$PRJ/src/lmps/thermal_cycle.lmp" 2>&1 | tee "$PRJ/logs/thermal_$name.log"

  mkdir -p "$out/sme/restarts"
  [[ -f "$out/sme/final.data" ]] || INPUT_DATA="unused" OUTPUT_DIR="$out/sme" RESTART_DIR="$out/sme/restarts" START_CMD="read_restart $out/thermal_cycle/restart_files/restart.stage2_end" POT_LIB="$LIB" POT_FILE="$POT" SEED=42 mpirun -n $MPI_RANKS "$LMP" -in "$PRJ/src/lmps/sme.lmp" 2>&1 | tee "$PRJ/logs/sme_$name.log"

  mkdir -p "$out/se/restarts"
  [[ -f "$out/se/final.data" ]] || INPUT_DATA="unused" OUTPUT_DIR="$out/se" RESTART_DIR="$out/se/restarts" THERMAL_RESTART="$out/thermal_cycle/restart_files/restart.stage4_end" POT_LIB="$LIB" POT_FILE="$POT" SEED=42 mpirun -n $MPI_RANKS "$LMP" -in "$PRJ/src/lmps/se.lmp" 2>&1 | tee "$PRJ/logs/se_$name.log"
done