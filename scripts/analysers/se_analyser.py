import glob
import re
import os
import sys
import numpy as np
import matplotlib
# FORCE HEADLESS MODE (Crucial for servers)
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, MaxNLocator
from scipy.signal import savgol_filter
from ovito.io import import_file
from ovito.modifiers import (
    PolyhedralTemplateMatchingModifier,
    ClusterAnalysisModifier,
    ExpressionSelectionModifier
)

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================
# Point this to your SE directory
INPUT_DIR = "/home/rimuru/workspace/data/raw_output/niti_10nm_Ni50_single_crystal/se"
OUTPUT_DIR = os.path.join(INPUT_DIR, "analysis_output")

if not os.path.exists(INPUT_DIR):
    print(f"ERROR: Input directory does not exist!\n{INPUT_DIR}")
    sys.exit(1)

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
    print(f"I: Created output directory: {OUTPUT_DIR}")

# ==============================================================================
# 2. PHYSICS PARAMETERS (Must match se.lmp)
# ==============================================================================
STRAIN_MAX = 0.10      # 10% Strain
DT         = 0.002     # Timestep (ps)
RMSD       = 0.15      # PTM Sensitivity
CORE_FRAC  = 0.8       # Definition of "Core" (r < 0.8*R)

# Protocol Timing (From se.lmp)
# load_time = 2000.0 -> 1,000,000 steps
# unload_time = 2000.0 -> 1,000,000 steps
# eq_steps = 1000
STEPS_LOAD   = 1000000
STEPS_UNLOAD = 1000000
EQ_OFFSET    = 1000

END_LOAD   = EQ_OFFSET + STEPS_LOAD
END_UNLOAD = END_LOAD + STEPS_UNLOAD

# ==============================================================================
# 3. ANALYSIS ENGINE
# ==============================================================================
def get_strain(step):
    """Maps timestep to Engineering Strain based on protocol."""
    if step < EQ_OFFSET:
        return 0.0
    elif step <= END_LOAD:
        # Loading: 0 -> 0.10
        return ((step - EQ_OFFSET) / STEPS_LOAD) * STRAIN_MAX
    elif step <= END_UNLOAD:
        # Unloading: 0.10 -> 0
        return STRAIN_MAX - ((step - END_LOAD) / STEPS_UNLOAD) * STRAIN_MAX
    else:
        return 0.0

print(f"I: Analyzing SE data in: {INPUT_DIR}")
# Look for dump.se.*.gz
files = sorted(glob.glob(os.path.join(INPUT_DIR, "dump.se.*.gz")), 
               key=lambda f: int(re.search(r'dump\.se\.(\d+)', f).group(1)))

if not files:
    print("ERROR: No dump files found (dump.se.*.gz).")
    sys.exit(1)

# --- MODIFIERS SETUP ---
print("I: Configuring OVITO Modifiers...")

# 1. PTM (Identify Phase)
ptm = PolyhedralTemplateMatchingModifier(rmsd_cutoff=RMSD, output_orientation=False)

# 2. Selection (Select NON-Austenite for Cluster Analysis)
# Structure Type 3 = BCC (Austenite in PTM for NiTi usually maps to BCC/B2)
select_mod = ExpressionSelectionModifier(expression = "StructureType != 3")

# 3. Cluster Analysis (Count Martensite Domains)
cluster_mod = ClusterAnalysisModifier(cutoff=3.5, sort_by_size=True, only_selected=True)

# Build Pipeline
pipeline = import_file(files[0])
pipeline.modifiers.append(ptm)
pipeline.modifiers.append(select_mod)
pipeline.modifiers.append(cluster_mod)

results = [] 

print(f"I: Processing {len(files)} frames (OVITO Engine)...")

for i, fname in enumerate(files):
    step = int(re.search(r'dump\.se\.(\d+)', fname).group(1))
    
    # Skip equilibration frames if you want cleaner plots, or keep them to show T stability
    # Here we keep everything but calculate strain carefully
    strain = get_strain(step)
    
    # Load Frame
    pipeline.source.load(fname)
    data = pipeline.compute()
    
    # --- A. GEOMETRY (Core Definition) ---
    pos = data.particles.position
    com = np.mean(pos, axis=0)
    r_sq = np.sum((pos - com)**2, axis=1)
    max_r2 = np.max(r_sq)
    radius = np.sqrt(max_r2) # Approximate Outer Radius
    
    is_core = r_sq < (CORE_FRAC * radius)**2
    is_surf = ~is_core
    
    # --- B. PHASES (Austenite = Type 3) ---
    structs = data.particles['Structure Type']
    n_total = len(structs)
    
    # Global
    f_aus = np.count_nonzero(structs == 3) / n_total
    
    # Core vs Surface Deconvolution
    if np.sum(is_core) > 0:
        f_aus_c = np.count_nonzero(structs[is_core] == 3) / np.sum(is_core)
    else:
        f_aus_c = 0.0
        
    if np.sum(is_surf) > 0:
        f_aus_s = np.count_nonzero(structs[is_surf] == 3) / np.sum(is_surf)
    else:
        f_aus_s = 0.0
    
    # --- C. STRESS (Virial -> GPa) ---
    # c_st is Pressure*Volume (bars*Ang^3) in 'metal' units
    # We want Sigma_ZZ for the CORE to avoid surface fluctuations
    stress_gpa = 0.0
    if 'c_st' in data.particles:
        st_tensor = data.particles['c_st']
        # Extract ZZ component (Index 2 in 0-based array: xx, yy, zz, yz, xz, xy)
        s_zz_accum = np.sum(st_tensor[is_core, 2]) 
        
        # Calculate Core Volume (Approximate as Sphere for Normalization)
        vol_core = (4.0/3.0) * np.pi * ((radius * CORE_FRAC)**3)
        
        # Convert: (bar*Ang^3) / Ang^3 = bar. 1 bar = 0.0001 GPa
        # Note: 'stress/atom' is negative of stress tensor in some conventions, 
        # but usually consistent. Positive = Tension, Negative = Compression.
        # SE is compression -> we expect Negative Stress.
        stress_gpa = (s_zz_accum / vol_core) * 0.0001

    # --- D. ENERGY ---
    pe_core = 0.0
    if 'c_pe_atom' in data.particles:
        pe_arr = data.particles['c_pe_atom']
        if np.sum(is_core) > 0:
            pe_core = np.mean(pe_arr[is_core])
            
    # --- E. MICROSTRUCTURE ---
    n_clusters = 0
    if 'ClusterAnalysis.cluster_count' in data.attributes:
        n_clusters = data.attributes['ClusterAnalysis.cluster_count']
        
    # --- F. GEOMETRIC BREATHING (Rg) ---
    rg = np.sqrt(np.mean(r_sq))

    # --- G. THERMAL CHECK ---
    # Assuming 'temp' is not in dump (unless computed), we skip or infer.
    # Actually, c_pe is a proxy for T if KE were included, but we'll use PE as thermodynamic proxy.
    
    results.append([strain, stress_gpa, f_aus, f_aus_c, f_aus_s, pe_core, n_clusters, rg, step])
    
    if i % 20 == 0: 
        print(f"   Step {step}: Strain={strain*100:.1f}% | Stress={stress_gpa:.2f} GPa | Aus={f_aus:.2f}")

# Convert & Sort by Step
arr = np.array(results)
arr = arr[arr[:,8].argsort()]

# Save raw data
np.savetxt(f"{OUTPUT_DIR}/se_extracted_data.txt", arr, 
           header="Strain Stress_Core_GPa f_Aus f_Aus_Core f_Aus_Surf PE_Core Clusters Rg Step")

# ==============================================================================
# 4. PLOTTING ENGINE (SE DASHBOARD STYLE)
# ==============================================================================
print("I: Generating SE Dashboard Plots...")

# --- A. STYLE CONFIGURATION ---
C_LOAD_RAW   = '#A6C9E2'   # Faint Blue
C_LOAD_LINE  = '#004C99'   # Strong Blue (Loading)
C_UNLOAD_RAW = '#F4A6A6'   # Faint Red
C_UNLOAD_LINE= '#D62828'   # Strong Red (Unloading)
C_FILL       = '#E0E0E0'   # Hysteresis Fill

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 16,
    'axes.linewidth': 2.5,
    'lines.linewidth': 3.0,
    'grid.alpha': 0.6,
})

# --- B. DATA SPLITTING ---
# Split into Loading and Unloading phases based on Step
load_mask = (arr[:, 8] <= END_LOAD) & (arr[:, 8] >= EQ_OFFSET)
unload_mask = (arr[:, 8] > END_LOAD)

d_load = arr[load_mask]
d_unload = arr[unload_mask]

# --- C. UNIVERSAL PLOTTING FUNCTION ---
def make_hysteresis_plot(ax, x_l, y_l, x_u, y_u, ylab, title, show_legend=False, invert_y=False):
    # 1. Hysteresis Fill
    if len(x_l) > 1 and len(x_u) > 1:
        # Interpolate to common grid for filling area
        common_x = np.linspace(0, STRAIN_MAX, 100)
        yl_interp = np.interp(common_x, x_l, y_l)
        yu_interp = np.interp(common_x, x_u[::-1], y_u[::-1]) # Unload is reversed in time
        ax.fill_between(common_x, yl_interp, yu_interp, color=C_FILL, alpha=0.3)

    # 2. Raw Data (Ghost)
    ax.plot(x_l, y_l, 'o', color=C_LOAD_LINE, markersize=2, alpha=0.1)
    ax.plot(x_u, y_u, 'o', color=C_UNLOAD_LINE, markersize=2, alpha=0.1)

    # 3. Smooth Trend Lines
    win = 11 if len(x_l) > 11 else 3
    yl_sm = savgol_filter(y_l, win, 2) if len(x_l) > win else y_l
    yu_sm = savgol_filter(y_u, win, 2) if len(x_u) > win else y_u
    
    ax.plot(x_l, yl_sm, '-', color=C_LOAD_LINE, label=r'Loading ($\varepsilon \uparrow$)', alpha=0.9)
    ax.plot(x_u, yu_sm, '-', color=C_UNLOAD_LINE, label=r'Unloading ($\varepsilon \downarrow$)', alpha=0.9)
    
    # 4. Formatting
    ax.set_xlabel("Compressive Strain", fontweight='bold')
    ax.set_ylabel(ylab, fontweight='bold')
    ax.set_title(title, fontweight='bold', fontsize=18, pad=12)
    
    ax.xaxis.set_major_locator(MultipleLocator(0.02)) # 2% steps
    ax.xaxis.set_minor_locator(MultipleLocator(0.005))
    
    ax.grid(True, which='major', linestyle='-', color='#888888', alpha=0.5)
    
    if invert_y:
        ax.invert_yaxis()
        
    if show_legend:
        ax.legend(fontsize=12, loc='best', frameon=True)

# --- D. INDIVIDUAL PLOTS ---
plots_config = [
    # (Filename, Col_Idx, YLabel, Title, InvertY?)
    ("1_Stress_Strain", 1, "Stress (GPa)", "Constitutive Response", True), # Compressive stress is negative
    ("2_Phase_Fraction", 2, "Austenite Fraction", "Phase Transformation", False),
    ("3_Core_Energy", 5, "PE Core (eV/atom)", "Energy Dissipation", True), # Lower Energy is better
    ("4_Cluster_Count", 6, "Martensite Clusters", "Nucleation Kinetics", False),
    ("5_Breathing", 7, r"Radius of Gyration ($\AA$)", "Geometric Breathing", False),
]

for fname, idx, ylab, title, inv_y in plots_config:
    fig = plt.figure(figsize=(10, 8))
    make_hysteresis_plot(plt.gca(), d_load[:,0], d_load[:,idx], d_unload[:,0], d_unload[:,idx], ylab, title, True, inv_y)
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/SE_{fname}.png", dpi=300)
    plt.close(fig)

# --- E. SPECIAL PLOT: CORE vs SURFACE LAG ---
fig, ax = plt.subplots(figsize=(10, 8))
# Smooth lines
l_core = savgol_filter(d_load[:,3], 11, 2)
l_surf = savgol_filter(d_load[:,4], 11, 2)
u_core = savgol_filter(d_unload[:,3], 11, 2)
u_surf = savgol_filter(d_unload[:,4], 11, 2)

ax.plot(d_load[:,0], l_core, '-', color='#228B22', lw=3, label='Core (Load)')
ax.plot(d_load[:,0], l_surf, '--', color='#800080', lw=3, label='Surface (Load)')
ax.plot(d_unload[:,0], u_core, '-', color='#228B22', alpha=0.4, lw=2)
ax.plot(d_unload[:,0], u_surf, '--', color='#800080', alpha=0.4, lw=2)

ax.set_xlabel("Compressive Strain", fontweight='bold')
ax.set_ylabel("Austenite Fraction", fontweight='bold')
ax.set_title("Nanoscale Size Effect: Core vs Surface", fontweight='bold', fontsize=18)
ax.legend(fontsize=14)
ax.grid(True)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/SE_Core_vs_Surface.png", dpi=300)
plt.close(fig)

# --- F. MASTER DASHBOARD ---
print("I: Compiling SE Master Dashboard...")
fig, axes = plt.subplots(2, 3, figsize=(28, 18))
fig.suptitle(f"NiTi Superelasticity Analysis Suite", fontsize=28, fontweight='bold', y=0.95)
ax_list = axes.flatten()

# 1. Stress-Strain (The Money Plot)
make_hysteresis_plot(ax_list[0], d_load[:,0], d_load[:,1], d_unload[:,0], d_unload[:,1], 
                     "Stress (GPa)", "Mechanical Response", True, True)

# 2. Phase Fraction
make_hysteresis_plot(ax_list[1], d_load[:,0], d_load[:,2], d_unload[:,0], d_unload[:,2], 
                     "Austenite Frac", "Phase Transformation", False, False)

# 3. Core vs Surface (Custom)
ax = ax_list[2]
ax.plot(d_load[:,0], l_core, '-', color='#228B22', lw=3, label='Core')
ax.plot(d_load[:,0], l_surf, '--', color='#800080', lw=3, label='Surface')
ax.set_title("Core vs Surface Lag", fontweight='bold', fontsize=18)
ax.set_xlabel("Strain")
ax.grid(True)
ax.legend()

# 4. Energy Dissipation
make_hysteresis_plot(ax_list[3], d_load[:,0], d_load[:,5], d_unload[:,0], d_unload[:,5], 
                     "PE (eV)", "Thermodynamics", False, True)

# 5. Nucleation (Clusters)
make_hysteresis_plot(ax_list[4], d_load[:,0], d_load[:,6], d_unload[:,0], d_unload[:,6], 
                     "Cluster Count", "Nucleation Density", False, False)

# 6. Combined (Overlay Stress + Phase)
ax6 = ax_list[5]
ax6.set_title("Constitutive Coupling", fontweight='bold', fontsize=18)
ax6.set_xlabel("Strain")
ax6_r = ax6.twinx()

ax6.plot(d_load[:,0], savgol_filter(d_load[:,1], 11, 2), 'b-', lw=3, label='Stress')
ax6_r.plot(d_load[:,0], savgol_filter(d_load[:,2], 11, 2), 'r--', lw=3, label='Phase')

ax6.set_ylabel("Stress (GPa)", color='b', fontweight='bold')
ax6_r.set_ylabel("Austenite Frac", color='r', fontweight='bold')
ax6.invert_yaxis() # Compressive stress down
ax6.grid(True)

plt.subplots_adjust(wspace=0.25, hspace=0.3)
plt.savefig(f"{OUTPUT_DIR}/SE_MASTER_DASHBOARD.png", dpi=300)

print(f"I: Analysis Complete. Data and Plots saved to: {OUTPUT_DIR}")