import glob
import re
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import savgol_filter
from ovito.io import import_file
from ovito.modifiers import (
    PolyhedralTemplateMatchingModifier,
    ClusterAnalysisModifier,
    AtomicStrainModifier
)

# ==============================================================================
# 1. USER CONFIGURATION (HARDCODE YOUR PATH HERE)
# ==============================================================================
# POINT THIS TO YOUR RAW DATA FOLDER
INPUT_DIR = "/home/rimuru/workspace/data/raw_output/niti_10nm_Ni50_single_crystal/thermal_cycle"

# AUTOMATIC OUTPUT PATH GENERATION
# Flips 'raw_output' to 'analysis_output' automatically
OUTPUT_DIR = INPUT_DIR.replace("raw_output", "analysis_output")

# Check Paths
if not os.path.exists(INPUT_DIR):
    print(f"ERROR: Input directory does not exist!\n{INPUT_DIR}")
    sys.exit(1)

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
    print(f"I: Created output directory: {OUTPUT_DIR}")

# ==============================================================================
# 2. JOURNAL STYLE SETTINGS (Acta Materialia / Scripta Mat)
# ==============================================================================
# Clean, professional, high-contrast
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['legend.frameon'] = False

# "Nature" Style Color Palette
C_COOL = '#0072B2'  # Blue
C_HEAT = '#D55E00'  # Red
C_CORE = '#009E73'  # Green
C_SURF = '#CC79A7'  # Purple

# ==============================================================================
# 3. PHYSICS PARAMETERS (MATCH LAMMPS)
# ==============================================================================
T_LOW  = 200.0
T_HIGH = 500.0
RATE   = 0.35
DT     = 0.002
RMSD   = 0.15   # Critical for NiTi PTM
CORE_FRAC = 0.8 # Core radius fraction

# Timing
DELTA_T    = T_HIGH - T_LOW
RUN_RAMP   = round(DELTA_T / RATE / DT)
RUN_HOLD   = 50000
END_STAGE0 = RUN_HOLD
END_STAGE1 = END_STAGE0 + RUN_RAMP
END_STAGE2 = END_STAGE1 + RUN_HOLD
END_STAGE3 = END_STAGE2 + RUN_RAMP

# ==============================================================================
# 4. ANALYSIS ENGINE
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

# --- MODIFIERS ---
ptm = PolyhedralTemplateMatchingModifier(rmsd_cutoff=RMSD, output_orientation=False)
ref_file = files[0] # Use first frame as reference for strain
strain_mod = AtomicStrainModifier(cutoff=3.5, reference=import_file(ref_file))
cluster = ClusterAnalysisModifier(cutoff=3.5, sort_by_size=True, only_selected=True)

results = [] 

print(f"I: Processing {len(files)} frames (OVITO Engine)...")

for i, fname in enumerate(files):
    step = int(re.search(r'dump\.cycle\.(\d+)', fname).group(1))
    temp = get_temp(step)
    
    pipeline = import_file(fname)
    pipeline.modifiers.append(ptm)
    pipeline.modifiers.append(strain_mod) 
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
    selection = data.particles.create_property("Selection")
    selection[...] = (structs != 3) # Select NON-Austenite
    pipeline.modifiers.append(cluster)
    data_cluster = pipeline.compute()
    n_clusters = data_cluster.attributes['ClusterAnalysis.cluster_count']
    
    # --- F. BREATHING ---
    rg = np.sqrt(np.mean(r_sq))
    
    # [Temp, f_aus, f_aus_c, f_aus_s, pe_c, shear, stress, clusters, rg, step]
    results.append([temp, f_aus, f_aus_c, f_aus_s, pe_core, shear_avg, stress_gpa, n_clusters, rg, step])
    
    if i % 50 == 0: print(f"   Step {step}: T={temp:.0f}K | Aus={f_aus:.2f}")

# Convert & Sort
arr = np.array(results)
arr = arr[arr[:,9].argsort()]

# Save raw data for Excel/Origin
np.savetxt(f"{OUTPUT_DIR}/extracted_data.txt", arr, header="Temp f_Aus f_Aus_Core f_Aus_Surf PE_Core Shear Stress Clusters Rg Step")

# ==============================================================================
# 5. PLOTTING ENGINE (JOURNAL READY)
# ==============================================================================
cool_mask = arr[:, 9] <= END_STAGE2
heat_mask = arr[:, 9] > END_STAGE2
cool = arr[cool_mask]
heat = arr[heat_mask]

def make_plot(ax, yc, yh, ylab, title, show_legend=False):
    ax.plot(cool[:,0], yc, '-', color=C_COOL, lw=2, label='Cooling ($A \\to M$)')
    ax.plot(heat[:,0], yh, '-', color=C_HEAT, lw=2, label='Heating ($M \\to A$)')
    ax.set_xlabel('Temperature (K)', fontweight='bold')
    ax.set_ylabel(ylab, fontweight='bold')
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.4)
    if show_legend: ax.legend(loc='best', fontsize=9)

print("I: Generating Figures...")

# --- 1. INDIVIDUAL GRAPHS ---
plots_config = [
    ("01_Hysteresis_Global", arr[:,1], "Austenite Fraction ($f_A$)", "Global Hysteresis"),
    ("02_Hysteresis_Core",   arr[:,2], "Core Austenite Fraction",     "Core Transformation"),
    ("03_Energy_DSC",        arr[:,4], "Potential Energy (eV/atom)",  "Thermodynamic Signature"),
    ("04_Shear_Strain",      arr[:,5], "Avg. Shear Strain",           "Lattice Distortion"),
    ("05_Internal_Stress",   arr[:,6], "Hydrostatic Stress (GPa)",    "Internal Conf. Stress"),
    ("06_Nucleation",        arr[:,7], "Cluster Count ($N_c$)",       "Nucleation Density"),
    ("07_Breathing_Rg",      arr[:,8], "Radius of Gyration ($\AA$)",  "Geometric Breathing")
]

for fname, col_data, ylab, title in plots_config:
    plt.figure(figsize=(6,5))
    yc = col_data[cool_mask]
    yh = col_data[heat_mask]
    make_plot(plt.gca(), yc, yh, ylab, title, show_legend=True)
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/{fname}.png", dpi=300)
    plt.close()

# --- 2. MASTER DASHBOARD (ALL SUBPLOTS) ---
fig, axes = plt.subplots(3, 3, figsize=(18, 15))
fig.suptitle(f"Thermal Characterization Suite", fontsize=16, fontweight='bold')

# Row 1
make_plot(axes[0,0], cool[:,1], heat[:,1], "$f_{Aus}$ (Global)", "(a) Global Hysteresis", True)
make_plot(axes[0,1], cool[:,2], heat[:,2], "$f_{Aus}$ (Core)",   "(b) Core Hysteresis")
make_plot(axes[0,2], cool[:,4], heat[:,4], "PE (eV/atom)",       "(c) Thermodynamics")

# Row 2
make_plot(axes[1,0], cool[:,5], heat[:,5], "Shear Strain",       "(d) Shear Order Parameter")
make_plot(axes[1,1], cool[:,6], heat[:,6], "Stress (GPa)",       "(e) Internal Stress")
make_plot(axes[1,2], cool[:,8], heat[:,8], "Rg ($\AA$)",         "(f) Geometric Breathing")

# Row 3
make_plot(axes[2,0], cool[:,7], heat[:,7], "Cluster Count",      "(g) Nucleation Density")

# (h) Core vs Surface Lag
axes[2,1].plot(cool[:,0], cool[:,2], '-', color=C_CORE, lw=2, label='Core')
axes[2,1].plot(cool[:,0], cool[:,3], '--', color=C_SURF, lw=2, label='Surface')
axes[2,1].set_xlabel('Temperature (K)', fontweight='bold')
axes[2,1].set_ylabel('$f_{Aus}$', fontweight='bold')
axes[2,1].set_title('(h) Surface vs Core Lag', fontsize=11, fontweight='bold')
axes[2,1].legend()
axes[2,1].grid(True, linestyle='--', alpha=0.4)

# (i) Kinetics
d_frac = np.gradient(cool[:,1], cool[:,0]) 
axes[2,2].plot(cool[:,0], np.abs(d_frac), '-', color=C_COOL, lw=1.5)
axes[2,2].set_xlabel('Temperature (K)', fontweight='bold')
axes[2,2].set_ylabel('|df/dT|', fontweight='bold')
axes[2,2].set_title('(i) Transformation Rate', fontsize=11, fontweight='bold')
axes[2,2].grid(True, linestyle='--', alpha=0.4)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(f"{OUTPUT_DIR}/MASTER_DASHBOARD.png", dpi=300)

print(f"I: Success! All graphs saved to: {OUTPUT_DIR}")