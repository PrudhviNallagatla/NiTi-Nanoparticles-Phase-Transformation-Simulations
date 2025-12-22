import glob
import re
import os
import sys
import numpy as np
import matplotlib
# FORCE HEADLESS MODE (Fixes QDBusError/Display crashes on Linux servers)
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
# Import MaxNLocator to fix the "Missing Grid/Y-values" bug
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, MaxNLocator
from scipy.signal import savgol_filter
from ovito.io import import_file
from ovito.modifiers import (
    PolyhedralTemplateMatchingModifier,
    ClusterAnalysisModifier,
    AtomicStrainModifier,
    ExpressionSelectionModifier
)

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================
INPUT_DIR = "/home/rimuru/workspace/data/raw_output/niti_10nm_Ni50_single_crystal/thermal_cycle"
OUTPUT_DIR = INPUT_DIR.replace("raw_output", "analysis_output")

if not os.path.exists(INPUT_DIR):
    print(f"ERROR: Input directory does not exist!\n{INPUT_DIR}")
    sys.exit(1)

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
    print(f"I: Created output directory: {OUTPUT_DIR}")

# ==============================================================================
# 2. PHYSICS PARAMETERS
# ==============================================================================
T_LOW  = 200.0
T_HIGH = 500.0
RATE   = 0.35
DT     = 0.002
RMSD   = 0.15
CORE_FRAC = 0.8

# Timing
DELTA_T    = T_HIGH - T_LOW
RUN_RAMP   = round(DELTA_T / RATE / DT)
RUN_HOLD   = 50000
END_STAGE0 = RUN_HOLD
END_STAGE1 = END_STAGE0 + RUN_RAMP
END_STAGE2 = END_STAGE1 + RUN_HOLD
END_STAGE3 = END_STAGE2 + RUN_RAMP

# ==============================================================================
# 3. ANALYSIS ENGINE
# ==============================================================================
def get_temp(step):
    if step <= END_STAGE0: return T_HIGH
    elif step <= END_STAGE1: return T_HIGH - ((step-END_STAGE0)/RUN_RAMP)*(T_HIGH-T_LOW)
    elif step <= END_STAGE2: return T_LOW
    elif step <= END_STAGE3: return T_LOW + ((step-END_STAGE2)/RUN_RAMP)*(T_HIGH-T_LOW)
    else: return T_HIGH

print(f"I: Analyzing data in: {INPUT_DIR}")
files = sorted(glob.glob(os.path.join(INPUT_DIR, "dump.cycle.*.gz")), 
               key=lambda f: int(re.search(r'dump\.cycle\.(\d+)', f).group(1)))

if not files:
    print("ERROR: No dump files found.")
    sys.exit(1)

# --- MODIFIERS SETUP ---
print("I: Configuring OVITO Modifiers...")

# 1. PTM
ptm = PolyhedralTemplateMatchingModifier(rmsd_cutoff=RMSD, output_orientation=False)

# 2. Strain (Requires Reference Pipeline Source)
ref_pipeline = import_file(files[0]) 
strain_mod = AtomicStrainModifier(cutoff=3.5)
strain_mod.reference = ref_pipeline.source 

# 3. Selection (Select NON-Austenite)
# Assuming Structure Type 3 is Austenite (B2)
select_mod = ExpressionSelectionModifier(expression = "StructureType != 3")

# 4. Cluster (Runs on selection)
cluster_mod = ClusterAnalysisModifier(cutoff=3.5, sort_by_size=True, only_selected=True)

# Build Pipeline
pipeline = import_file(files[0])
pipeline.modifiers.append(ptm)
pipeline.modifiers.append(strain_mod)
pipeline.modifiers.append(select_mod)
pipeline.modifiers.append(cluster_mod)

results = [] 

print(f"I: Processing {len(files)} frames (OVITO Engine)...")

for i, fname in enumerate(files):
    step = int(re.search(r'dump\.cycle\.(\d+)', fname).group(1))
    temp = get_temp(step)
    
    # Load Frame
    pipeline.source.load(fname)
    data = pipeline.compute()
    
    # --- A. GEOMETRY ---
    pos = data.particles.position
    com = np.mean(pos, axis=0)
    r_sq = np.sum((pos - com)**2, axis=1)
    max_r2 = np.max(r_sq)
    radius = np.sqrt(max_r2)
    is_core = r_sq < (CORE_FRAC * radius)**2
    is_surf = ~is_core
    
    # --- B. PHASES (Austenite = Type 3) ---
    structs = data.particles['Structure Type']
    n_total = len(structs)
    f_aus = np.count_nonzero(structs == 3) / n_total
    f_aus_c = np.count_nonzero(structs[is_core] == 3) / np.sum(is_core) if np.sum(is_core) > 0 else 0
    f_aus_s = np.count_nonzero(structs[is_surf] == 3) / np.sum(is_surf) if np.sum(is_surf) > 0 else 0
    
    # --- C. ENERGY ---
    if 'c_pe_all' in data.particles:
        pe_arr = data.particles['c_pe_all']
        pe_core = np.mean(pe_arr[is_core]) if np.sum(is_core) > 0 else 0
    else: pe_core = 0.0
        
    # --- D. MECHANICS ---
    shear = data.particles['Shear Strain']
    shear_avg = np.mean(shear)
    
    stress_gpa = 0.0
    if 'c_st' in data.particles:
        st_tensor = data.particles['c_st']
        trace_sum = np.sum((st_tensor[:,0] + st_tensor[:,1] + st_tensor[:,2])[is_core])
        vol_core = (4/3) * np.pi * (radius * CORE_FRAC)**3
        stress_gpa = (trace_sum / 3.0) / vol_core / 10000.0
        
    # --- E. MICROSTRUCTURE ---
    if 'ClusterAnalysis.cluster_count' in data.attributes:
        n_clusters = data.attributes['ClusterAnalysis.cluster_count']
    else:
        n_clusters = 0
    
    # --- F. BREATHING ---
    rg = np.sqrt(np.mean(r_sq))
    
    results.append([temp, f_aus, f_aus_c, f_aus_s, pe_core, shear_avg, stress_gpa, n_clusters, rg, step])
    
    if i % 50 == 0: print(f"   Step {step}: T={temp:.0f}K | Aus={f_aus:.2f}")

# Convert & Sort
arr = np.array(results)
arr = arr[arr[:,9].argsort()]

# Save raw data
np.savetxt(f"{OUTPUT_DIR}/extracted_data.txt", arr, header="Temp f_Aus f_Aus_Core f_Aus_Surf PE_Core Shear Stress Clusters Rg Step")

# ==============================================================================
# 4. PLOTTING ENGINE (FULL SIZE & CORRECT NAMES)
# ==============================================================================
print("I: Generating High-Precision Plots...")

# --- A. STYLE CONFIGURATION ---
C_COOL_RAW  = '#A6C9E2'  # Faint Blue (Raw Data)
C_COOL_LINE = '#004C99'  # Strong Blue (Trend)
C_HEAT_RAW  = '#F4A6A6'  # Faint Red (Raw Data)
C_HEAT_LINE = '#D62828'  # Strong Red (Trend)
C_FILL      = '#E0E0E0'  # Grey (Hysteresis Area)

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 16,              # Bigger font
    'axes.linewidth': 2.5,
    'lines.linewidth': 3.0,
    'xtick.major.size': 10,
    'ytick.major.size': 10,
    'grid.alpha': 0.6,
})

# --- B. DATA PREPARATION ---
cool_mask = (arr[:, 9] <= END_STAGE2) & (arr[:, 0] <= T_HIGH)
heat_mask = (arr[:, 9] > END_STAGE2) & (arr[:, 0] <= T_HIGH)

cool = arr[cool_mask]
heat = arr[heat_mask]
cool = cool[cool[:, 0].argsort()[::-1]] 
heat = heat[heat[:, 0].argsort()]       

# --- C. UNIVERSAL PLOTTING FUNCTION ---
def make_detailed_plot(ax, x_c, y_c, x_h, y_h, ylab, title, show_legend=False):
    # 1. Hysteresis Fill
    if len(x_c) > 1 and len(x_h) > 1:
        x_union = np.sort(np.unique(np.concatenate((x_c, x_h))))
        yc_interp = np.interp(x_union, x_c[::-1], y_c[::-1]) 
        yh_interp = np.interp(x_union, x_h, y_h)
        ax.fill_between(x_union, yc_interp, yh_interp, color=C_FILL, alpha=0.3)

    # 2. Ghost Plotting (Raw Data Dots)
    ax.plot(x_c, y_c, 'o', color=C_COOL_LINE, markersize=3, alpha=0.2, markeredgewidth=0)
    ax.plot(x_h, y_h, 'o', color=C_HEAT_LINE, markersize=3, alpha=0.2, markeredgewidth=0)

    # 3. Smooth Trend Lines
    win = 9 if len(x_c) > 9 else 3
    yc_sm = savgol_filter(y_c, win, 3) if len(x_c) > win else y_c
    yh_sm = savgol_filter(y_h, win, 3) if len(x_h) > win else y_h
    
    ax.plot(x_c, yc_sm, '-', color=C_COOL_LINE, label=r'Cooling ($A \to M$)', alpha=0.9)
    ax.plot(x_h, yh_sm, '-', color=C_HEAT_LINE, label=r'Heating ($M \to A$)', alpha=0.9)
    
    # 4. Pinpoint Grids
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    ax.grid(True, which='major', linestyle='-', color='#888888', alpha=0.5)
    ax.grid(True, which='minor', linestyle=':', color='#AAAAAA', alpha=0.4)

    ax.set_xlabel("Temperature (K)", fontweight='bold')
    ax.set_ylabel(ylab, fontweight='bold')
    ax.set_title(title, fontweight='bold', fontsize=20, pad=15)
    
    if show_legend: 
        ax.legend(fontsize=14, loc='best', frameon=True, framealpha=0.95)

# --- D. GENERATE INDIVIDUAL PLOTS ---
# Format: (Filename, Column_Index, Y_Label, Title)
plots_config = [
    ("Global_Hysteresis_Loop", 1, r"Austenite Fraction ($f_{Aus}$)", "Global Hysteresis Loop"),
    ("Core_Transformation",    2, r"Core Austenite ($f_{Core}$)", "Core Transformation"),
    ("Potential_Energy",       4, r"PE (eV/atom)", "Thermodynamic Response"), 
    ("Shear_Strain",           5, r"Avg. Shear Strain ($\eta_{VM}$)", "Lattice Distortion"),
    ("Internal_Stress",        6, r"Stress (GPa)", "Internal Stress Evolution"),
    ("Cluster_Analysis",       7, r"Cluster Count", "Nucleation Density"),
    ("Geometric_Breathing",    8, r"Radius of Gyration ($\AA$)", "Geometric Breathing"),
]

for fname, idx, ylab, title in plots_config:
    # Large, spacious individual plots
    fig = plt.figure(figsize=(12, 10)) 
    make_detailed_plot(plt.gca(), cool[:,0], cool[:,idx], heat[:,0], heat[:,idx], ylab, title, show_legend=True)
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/{fname}.png", dpi=300)
    plt.close(fig)

# --- E. MASTER DASHBOARD (FULL SIZE) ---
print("I: Compiling Master Dashboard (Full Size)...")
# Massive size (30 inches wide) to avoid "compacting"
fig, axes = plt.subplots(2, 3, figsize=(30, 20))
fig.suptitle(f"NiTi Thermal Analysis Suite", fontsize=30, fontweight='bold', y=0.96)
ax_list = axes.flatten()

# 1. Main Hysteresis
make_detailed_plot(ax_list[0], cool[:,0], cool[:,1], heat[:,0], heat[:,1], r"$f_{Aus}$", "Global Phase Fraction", True)

# 2. Core vs Surface
ax = ax_list[1]
ax.plot(cool[:,0], cool[:,2], '-', color='#228B22', label='Core', lw=4)
ax.plot(cool[:,0], cool[:,3], '--', color='#800080', label='Surface', lw=4)
ax.set_title("Core vs Surface Lag", fontweight='bold', fontsize=20)
ax.set_xlabel("Temperature (K)", fontweight='bold')
ax.set_ylabel("Austenite Fraction", fontweight='bold')
ax.xaxis.set_major_locator(MultipleLocator(20))
ax.grid(True, alpha=0.4)
ax.legend(fontsize=14)
ax.invert_xaxis()

# 3. Energy
make_detailed_plot(ax_list[2], cool[:,0], cool[:,4], heat[:,0], heat[:,4], "PE (eV)", "Thermodynamics")

# 4. Stress
make_detailed_plot(ax_list[3], cool[:,0], cool[:,6], heat[:,0], heat[:,6], "Stress (GPa)", "Internal Stress")

# 5. Breathing (Rg) - Replaced Clusters with Breathing per your interest
make_detailed_plot(ax_list[4], cool[:,0], cool[:,8], heat[:,0], heat[:,8], r"$R_g$ ($\AA$)", "Geometric Breathing")

# 6. Kinetics
ax = ax_list[5]
valid_mask = np.insert(np.diff(cool[:,0]) != 0, 0, True) 
T_clean = cool[valid_mask, 0]
f_clean = cool[valid_mask, 1]

if len(T_clean) > 10:
    grad = np.gradient(f_clean, T_clean)
    abs_grad = savgol_filter(np.abs(grad), 11, 2) if len(grad) > 11 else np.abs(grad)
    ax.plot(T_clean, abs_grad, color='#333333', lw=3)
    ax.fill_between(T_clean, 0, abs_grad, color='gray', alpha=0.2)
    ax.set_title("Transformation Rate |df/dT|", fontweight='bold', fontsize=20)
else:
    ax.text(0.5, 0.5, "Insufficient Data", ha='center', fontsize=18)

ax.set_xlabel("Temperature (K)", fontweight='bold')
ax.grid(True, alpha=0.4)
ax.invert_xaxis()

plt.subplots_adjust(left=0.05, right=0.95, top=0.92, bottom=0.05, wspace=0.25, hspace=0.25)
plt.savefig(f"{OUTPUT_DIR}/MASTER_DASHBOARD_FULL.png", dpi=300)

# ==============================================================================
# 5. TRANSFORMATION CALCULATOR
# ==============================================================================
print("\nI: Calculating Critical Temperatures...")

f_min, f_max = np.min(arr[:, 1]), np.max(arr[:, 1])
f_range = f_max - f_min

def get_temp_at_frac(T_arr, f_arr, threshold, mode='falling'):
    target = f_min + (threshold * f_range)
    f_sm = savgol_filter(f_arr, 11, 3) if len(f_arr) > 11 else f_arr
    if mode == 'falling': 
        idx = np.where(f_sm <= target)[0]
        return T_arr[idx[0]] if len(idx) > 0 else None
    else: 
        idx = np.where(f_sm >= target)[0]
        return T_arr[idx[0]] if len(idx) > 0 else None

As = get_temp_at_frac(heat[:,0], heat[:,1], 0.05, 'rising')
Af = get_temp_at_frac(heat[:,0], heat[:,1], 0.95, 'rising')
Ms = get_temp_at_frac(cool[:,0], cool[:,1], 0.95, 'falling')
Mf = get_temp_at_frac(cool[:,0], cool[:,1], 0.05, 'falling')

print("-" * 40)
print(f"PHASE TRANSFORMATION REPORT")
print("-" * 40)
print(f"As: {As:.1f} K" if As else "As: Not Detected")
print(f"Af: {Af:.1f} K" if Af else "Af: Not Detected")
print(f"Ms: {Ms:.1f} K" if Ms else "Ms: Not Detected")
print(f"Mf: {Mf:.1f} K" if Mf else "Mf: Not Detected")
print("-" * 40)

# Verification Plot
fig = plt.figure(figsize=(14, 12))
ax = plt.gca()
yc_sm = savgol_filter(cool[:,1], 11, 3) if len(cool) > 11 else cool[:,1]
yh_sm = savgol_filter(heat[:,1], 11, 3) if len(heat) > 11 else heat[:,1]

ax.plot(cool[:,0], cool[:,1], color=C_COOL_LINE, alpha=0.2, lw=1)
ax.plot(heat[:,0], heat[:,1], color=C_HEAT_LINE, alpha=0.2, lw=1)
ax.plot(cool[:,0], yc_sm, color=C_COOL_LINE, lw=4, label='Cooling')
ax.plot(heat[:,0], yh_sm, color=C_HEAT_LINE, lw=4, label='Heating')

if Ms: ax.plot(Ms, f_min + 0.95*f_range, 'o', color=C_COOL_LINE, ms=12, mec='k', label=f'Ms ({Ms:.0f}K)')
if Mf: ax.plot(Mf, f_min + 0.05*f_range, 'X', color=C_COOL_LINE, ms=12, mec='k', label=f'Mf ({Mf:.0f}K)')
if As: ax.plot(As, f_min + 0.05*f_range, 'o', color=C_HEAT_LINE, ms=12, mec='k', label=f'As ({As:.0f}K)')
if Af: ax.plot(Af, f_min + 0.95*f_range, 'X', color=C_HEAT_LINE, ms=12, mec='k', label=f'Af ({Af:.0f}K)')

ax.set_title("Transformation Temperatures Verification", fontweight='bold', fontsize=24)
ax.set_xlabel("Temperature (K)", fontweight='bold')
ax.set_ylabel("Austenite Fraction", fontweight='bold')
ax.legend(fontsize=16)
ax.grid(True, alpha=0.5)

plt.savefig(f"{OUTPUT_DIR}/Transformation_Verification.png", dpi=300)
print(f"I: Saved plots to {OUTPUT_DIR}")