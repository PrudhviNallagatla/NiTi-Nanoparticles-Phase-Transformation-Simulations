import glob
import re
import os
import sys
import numpy as np
import matplotlib
import math

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d
from ovito.io import import_file
from ovito.modifiers import (
    PolyhedralTemplateMatchingModifier,
    ClusterAnalysisModifier,
    ExpressionSelectionModifier,
)

INPUT_DIR = "/home/rimuru/workspace/data/raw_output/niti_10nm_Ni50_single_crystal/se"
OUTPUT_DIR = INPUT_DIR.replace("raw_output", "analysis_output")

if not os.path.exists(INPUT_DIR):
    print(f"ERROR: Input directory does not exist!\n{INPUT_DIR}")
    sys.exit(1)

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
    print(f"I: Created output directory: {OUTPUT_DIR}")

STRAIN_MAX = 0.10
DT = 0.002
RMSD = 0.15
CORE_FRAC = 0.8

STEPS_LOAD = 1000000
STEPS_UNLOAD = 1000000
EQ_OFFSET = 1000

END_LOAD = EQ_OFFSET + STEPS_LOAD
END_UNLOAD = END_LOAD + STEPS_UNLOAD


def get_strain(step):
    """Maps timestep to Engineering Strain (compressive = positive)."""
    if step < EQ_OFFSET:
        return 0.0
    elif step <= END_LOAD:
        return ((step - EQ_OFFSET) / STEPS_LOAD) * STRAIN_MAX
    elif step <= END_UNLOAD:
        return STRAIN_MAX * (1.0 - (step - END_LOAD) / STEPS_UNLOAD)
    else:
        return 0.0


def _get_particle_prop_array(data, name):
    """Return particle property as a NumPy array if it exists, else None."""
    try:
        if name in data.particles:
            return np.asarray(data.particles[name])
    except Exception:
        pass
    return None


def _find_stress_zz_key(data):
    """Find the ZZ stress component key."""
    keys = list(getattr(data.particles, "keys", lambda: [])())
    if not keys:
        return None

    if "c_st[3]" in data.particles:
        return "c_st[3]"

    preferred = [
        "c_st_3",
        "c_st3",
        "stzz",
        "szz",
        "Stress.ZZ",
        "Stress ZZ",
        "stress_zz",
        "stresszz",
    ]
    for k in preferred:
        if k in data.particles:
            return k

    pat = re.compile(r"^c_st\[\s*3\s*\]$", re.IGNORECASE)
    for k in keys:
        if pat.match(str(k)):
            return k

    for k in keys:
        sk = str(k).lower()
        if "c_st" in sk and ("[3]" in sk or " 3" in sk or sk.endswith("3")):
            return k

    return None


def compute_sigma_zz_gpa(data, volume_ang3):
    """Computes global sigma_zz (GPa) from LAMMPS stress/atom dump fields."""
    if volume_ang3 <= 0:
        return 0.0

    zz_key = _find_stress_zz_key(data)
    if zz_key is not None:
        stzz = _get_particle_prop_array(data, zz_key)
        if stzz is not None:
            sigma_bar = -float(np.sum(stzz)) / (volume_ang3 + 1e-12)
            return sigma_bar * 1.0e-4

    st = _get_particle_prop_array(data, "c_st")
    if st is not None and st.ndim == 2 and st.shape[1] >= 3:
        sigma_bar = -float(np.sum(st[:, 2])) / (volume_ang3 + 1e-12)
        return sigma_bar * 1.0e-4

    return 0.0


print(f"I: Analyzing SE data in: {INPUT_DIR}")
files = sorted(
    glob.glob(os.path.join(INPUT_DIR, "dump.se.*.gz")),
    key=lambda f: int(re.search(r"dump\.se\.(\d+)", f).group(1)),
)

if not files:
    print("ERROR: No dump files found (dump.se.*.gz).")
    sys.exit(1)

print("I: Configuring OVITO Modifiers...")

ptm = PolyhedralTemplateMatchingModifier(rmsd_cutoff=RMSD, output_orientation=False)
select_mod = ExpressionSelectionModifier(expression="StructureType != 3")
cluster_mod = ClusterAnalysisModifier(cutoff=3.5, sort_by_size=True, only_selected=True)

pipeline = import_file(files[0])
pipeline.modifiers.append(ptm)
pipeline.modifiers.append(select_mod)
pipeline.modifiers.append(cluster_mod)

results = []

print(f"I: Processing {len(files)} frames (OVITO Engine)...")

_warned_stress_missing = False
_printed_props = False

for i, fname in enumerate(files):
    step = int(re.search(r"dump\.se\.(\d+)", fname).group(1))
    strain = get_strain(step)

    pipeline.source.load(fname)
    data = pipeline.compute()

    if not _printed_props:
        keys = list(data.particles.keys())
        stress_keys = [
            k for k in keys if "c_st" in str(k).lower() or "stress" in str(k).lower()
        ]
        print("I: OVITO particle stress-related keys:", stress_keys[:50])
        zz_key = _find_stress_zz_key(data)
        print("I: Using ZZ stress key:", zz_key)
        _printed_props = True

    pos = np.asarray(data.particles.position)
    com = np.mean(pos, axis=0)
    r_sq = np.sum((pos - com) ** 2, axis=1)
    max_r2 = float(np.max(r_sq))
    radius = math.sqrt(max_r2)

    is_core = r_sq < (CORE_FRAC * radius) ** 2
    is_surf = ~is_core

    vol_sphere = (4.0 / 3.0) * math.pi * (radius**3)

    structs = np.asarray(data.particles["Structure Type"])
    n_total = len(structs)

    f_aus = np.count_nonzero(structs == 3) / (n_total + 1e-12)
    f_aus_c = np.count_nonzero(structs[is_core] == 3) / (np.sum(is_core) + 1e-12)
    f_aus_s = np.count_nonzero(structs[is_surf] == 3) / (np.sum(is_surf) + 1e-12)

    stress_gpa = compute_sigma_zz_gpa(data, vol_sphere)
    if stress_gpa == 0.0 and not _warned_stress_missing:
        zz_key = _find_stress_zz_key(data)
        has_any = (zz_key is not None) or ("c_st" in data.particles)
        if not has_any:
            print("W: No usable stress fields found. Stress will remain 0.")
            _warned_stress_missing = True

    pe_core = 0.0
    if "c_pe_atom" in data.particles:
        pe_arr = np.asarray(data.particles["c_pe_atom"])
        if np.sum(is_core) > 0:
            pe_core = float(np.mean(pe_arr[is_core]))

    n_clusters = 0
    if "ClusterAnalysis.cluster_count" in data.attributes:
        n_clusters = int(data.attributes["ClusterAnalysis.cluster_count"])

    rg = float(np.sqrt(np.mean(r_sq)))

    results.append(
        [strain, stress_gpa, f_aus, f_aus_c, f_aus_s, pe_core, n_clusters, rg, step]
    )

    if i % 20 == 0:
        print(
            f"   Step {step}: Strain={strain*100:.1f}% | Stress={stress_gpa:.2f} GPa | Aus={f_aus:.2f}"
        )

arr = np.array(results)
arr = arr[arr[:, 8].argsort()]

print(f"\nI: Data summary:")
print(f"   Total frames: {len(arr)}")
print(f"   Strain range: {arr[:,0].min():.4f} to {arr[:,0].max():.4f}")
print(f"   Stress range: {arr[:,1].min():.4f} to {arr[:,1].max():.4f} GPa")
print(f"   Austenite fraction: {arr[:,2].min():.4f} to {arr[:,2].max():.4f}")
print(f"   Step range: {int(arr[:,8].min())} to {int(arr[:,8].max())}")

np.savetxt(
    f"{OUTPUT_DIR}/se_extracted_data.txt",
    arr,
    header="Strain Stress_GPa f_Aus f_Aus_Core f_Aus_Surf PE_Core Clusters Rg Step",
)

C_LOAD = "#1f77b4"
C_UNLOAD = "#d62728"
C_CORE = "#2ca02c"
C_SURF = "#9467bd"
C_GRAY = "#444444"

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 11,
        "axes.labelsize": 13,
        "axes.titlesize": 14,
        "axes.titleweight": "bold",
        "axes.linewidth": 1.3,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "legend.fontsize": 10,
        "legend.frameon": True,
        "legend.framealpha": 0.9,
        "lines.linewidth": 1.8,  # REDUCED from 2.5
    }
)


def style_ax(ax):
    ax.grid(True, which="major", alpha=0.2, linewidth=0.9)
    ax.grid(True, which="minor", alpha=0.12, linewidth=0.6)
    ax.minorticks_on()
    ax.spines["top"].set_linewidth(1.3)
    ax.spines["right"].set_linewidth(1.3)
    ax.spines["left"].set_linewidth(1.3)
    ax.spines["bottom"].set_linewidth(1.3)


def smooth(y, window=7, poly=2):
    """Conservative smoothing."""
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 5:
        return y

    window = int(window)
    if window % 2 == 0:
        window += 1

    window = min(window, n if n % 2 == 1 else n - 1)
    window = max(window, poly + 2)

    if window > n or window < poly + 2:
        return y

    try:
        return savgol_filter(y, window, poly, mode="interp")
    except Exception as e:
        print(f"W: Smoothing failed: {e}")
        return y


# Split data
load_mask = (arr[:, 8] >= EQ_OFFSET) & (arr[:, 8] <= END_LOAD)
unload_mask = (arr[:, 8] > END_LOAD) & (arr[:, 8] <= END_UNLOAD)

d_load = arr[load_mask]
d_unload = arr[unload_mask]

print(f"I: Load phase: {len(d_load)} frames")
print(f"I: Unload phase: {len(d_unload)} frames")

if len(d_load) < 10 or len(d_unload) < 10:
    print("ERROR: Insufficient data points.")
    sys.exit(1)

d_load = d_load[d_load[:, 0].argsort()]
d_unload = d_unload[d_unload[:, 0].argsort()]

# Extract variables
eps_load = d_load[:, 0]
eps_unload = d_unload[:, 0]
sig_load = np.abs(d_load[:, 1])
sig_unload = np.abs(d_unload[:, 1])
phi_load = d_load[:, 2]
phi_unload = d_unload[:, 2]
phi_core_load = d_load[:, 3]
phi_surf_load = d_load[:, 4]
phi_core_unload = d_unload[:, 3]
phi_surf_unload = d_unload[:, 4]
pe_load = d_load[:, 5]
pe_unload = d_unload[:, 5]
clusters_load = d_load[:, 6]
clusters_unload = d_unload[:, 6]
rg_load = d_load[:, 7]
rg_unload = d_unload[:, 7]
step_load = d_load[:, 8]
step_unload = d_unload[:, 8]

pdf_path = os.path.join(OUTPUT_DIR, "Journal_SE_Complete_Analysis.pdf")
print(f"\nI: Generating comprehensive PDF: {pdf_path}")

with PdfPages(pdf_path) as pdf:

    # GRAPH 1: Mechanical Response (The Money Plot)
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    style_ax(ax)

    ax.plot(
        eps_load * 100,
        smooth(sig_load, 7),
        color=C_LOAD,
        linewidth=1.8,
        label="Loading",
        zorder=3,
    )
    ax.plot(
        eps_unload * 100,
        smooth(sig_unload, 7),
        color=C_UNLOAD,
        linewidth=1.8,
        label="Unloading",
        zorder=3,
    )

    # Fill hysteresis
    common_eps = np.linspace(0, STRAIN_MAX, 200)
    sig_l_interp = np.interp(common_eps, eps_load, smooth(sig_load, 7))
    sig_u_interp = np.interp(common_eps, eps_unload, smooth(sig_unload, 7))
    ax.fill_between(
        common_eps * 100,
        sig_l_interp,
        sig_u_interp,
        color="gray",
        alpha=0.15,
        label="Hysteresis",
    )

    ax.set_xlabel("Compressive Strain, ε (%)", fontsize=14, fontweight="bold")
    ax.set_ylabel(
        "Compressive Stress, |σ$_{zz}$| (GPa)", fontsize=14, fontweight="bold"
    )
    ax.set_title(
        "Mechanical Response: Stress-Strain Hysteresis Loop",
        fontsize=16,
        fontweight="bold",
        pad=15,
    )
    ax.legend(loc="best", fontsize=12)
    ax.set_xlim(0, STRAIN_MAX * 100)
    plt.tight_layout(pad=2.0)
    pdf.savefig(fig, dpi=300)
    plt.close(fig)

    # GRAPH 2: Tangent Modulus (Stiffness Degradation)
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    style_ax(ax)

    # Calculate tangent modulus
    E_load = np.gradient(smooth(sig_load, 9), eps_load + 1e-12)
    E_unload = np.gradient(smooth(sig_unload, 9), eps_unload + 1e-12)
    E_load = gaussian_filter1d(E_load, sigma=1.5)
    E_unload = gaussian_filter1d(E_unload, sigma=1.5)

    ax.plot(eps_load * 100, E_load, color=C_LOAD, linewidth=1.8, label="Loading")
    ax.plot(
        eps_unload * 100, E_unload, color=C_UNLOAD, linewidth=1.8, label="Unloading"
    )
    ax.axhline(0, color=C_GRAY, linewidth=1.5, linestyle="--", alpha=0.7)

    ax.set_xlabel("Compressive Strain, ε (%)", fontsize=14, fontweight="bold")
    ax.set_ylabel("Tangent Stiffness, dσ/dε (GPa)", fontsize=14, fontweight="bold")
    ax.set_title(
        "Stiffness Degradation: Tangent Modulus Evolution",
        fontsize=16,
        fontweight="bold",
        pad=15,
    )
    ax.legend(loc="best", fontsize=12)
    ax.set_xlim(0, STRAIN_MAX * 100)
    plt.tight_layout(pad=2.0)
    pdf.savefig(fig, dpi=300)
    plt.close(fig)

    # GRAPH 3: Energy Dissipation (Thermodynamics)
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    style_ax(ax)

    ax.plot(
        eps_load * 100, smooth(pe_load, 7), color=C_LOAD, linewidth=1.8, label="Loading"
    )
    ax.plot(
        eps_unload * 100,
        smooth(pe_unload, 7),
        color=C_UNLOAD,
        linewidth=1.8,
        label="Unloading",
    )

    # Fill energy hysteresis
    pe_l_interp = np.interp(common_eps, eps_load, smooth(pe_load, 7))
    pe_u_interp = np.interp(common_eps, eps_unload, smooth(pe_unload, 7))
    ax.fill_between(
        common_eps * 100,
        pe_l_interp,
        pe_u_interp,
        color="orange",
        alpha=0.2,
        label="Energy Dissipation",
    )

    ax.set_xlabel("Compressive Strain, ε (%)", fontsize=14, fontweight="bold")
    ax.set_ylabel("Potential Energy per Atom (eV)", fontsize=14, fontweight="bold")
    ax.set_title(
        "Energy Dissipation: Potential Energy Hysteresis",
        fontsize=16,
        fontweight="bold",
        pad=15,
    )
    ax.legend(loc="best", fontsize=12)
    ax.set_xlim(0, STRAIN_MAX * 100)
    plt.tight_layout(pad=2.0)
    pdf.savefig(fig, dpi=300)
    plt.close(fig)

    # GRAPH 4: Phase Transformation (The Cause) - IMPROVED
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    style_ax(ax)

    # Less smoothing for clarity
    ax.plot(
        eps_load * 100,
        smooth(phi_load, 5),
        color=C_LOAD,
        linewidth=1.8,
        label="Loading",
        marker="o",
        markersize=2.5,
        markevery=12,
    )
    ax.plot(
        eps_unload * 100,
        smooth(phi_unload, 5),
        color=C_UNLOAD,
        linewidth=1.8,
        label="Unloading",
        marker="s",
        markersize=2.5,
        markevery=12,
    )

    ax.set_xlabel("Compressive Strain, ε (%)", fontsize=14, fontweight="bold")
    ax.set_ylabel("Austenite Fraction, φ$_A$", fontsize=14, fontweight="bold")
    ax.set_title(
        "Phase Transformation: Austenite Fraction Evolution",
        fontsize=16,
        fontweight="bold",
        pad=15,
    )
    ax.legend(loc="best", fontsize=12)
    ax.set_xlim(0, STRAIN_MAX * 100)
    ax.set_ylim(0.2, 0.8)  # Changed from (0, 0.8)
    plt.tight_layout(pad=2.0)
    pdf.savefig(fig, dpi=300)
    plt.close(fig)

    # GRAPH 5: Transformation Kinetics (The Rate)
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    style_ax(ax)

    dphi_de_load = np.gradient(smooth(phi_load, 11), eps_load + 1e-12)
    dphi_de_unload = np.gradient(smooth(phi_unload, 11), eps_unload + 1e-12)
    dphi_de_load = gaussian_filter1d(dphi_de_load, sigma=1.5)
    dphi_de_unload = gaussian_filter1d(dphi_de_unload, sigma=1.5)

    ax.plot(
        eps_load * 100, dphi_de_load * 100, color=C_LOAD, linewidth=1.8, label="Loading"
    )
    ax.plot(
        eps_unload * 100,
        dphi_de_unload * 100,
        color=C_UNLOAD,
        linewidth=1.8,
        label="Unloading",
    )
    ax.axhline(0, color=C_GRAY, linewidth=1.5, linestyle="--", alpha=0.7)

    ax.set_xlabel("Compressive Strain, ε (%)", fontsize=14, fontweight="bold")
    ax.set_ylabel(
        "Transformation Rate, dφ$_A$/dε (%/strain)", fontsize=14, fontweight="bold"
    )
    ax.set_title(
        "Transformation Kinetics: Phase Evolution Rate",
        fontsize=16,
        fontweight="bold",
        pad=15,
    )
    ax.legend(loc="best", fontsize=12)
    ax.set_xlim(0, STRAIN_MAX * 100)
    plt.tight_layout(pad=2.0)
    pdf.savefig(fig, dpi=300)
    plt.close(fig)

    # GRAPH 6: Nucleation Density (Cluster Analysis)
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    style_ax(ax)

    ax.plot(
        eps_load * 100,
        smooth(clusters_load, 5),
        color=C_LOAD,
        linewidth=1.8,
        label="Loading",
        marker="o",
        markersize=3,
        markevery=10,
    )
    ax.plot(
        eps_unload * 100,
        smooth(clusters_unload, 5),
        color=C_UNLOAD,
        linewidth=1.8,
        label="Unloading",
        marker="s",
        markersize=3,
        markevery=10,
    )

    ax.set_xlabel("Compressive Strain, ε (%)", fontsize=14, fontweight="bold")
    ax.set_ylabel("Number of Martensite Clusters", fontsize=14, fontweight="bold")
    ax.set_title(
        "Nucleation Density: Martensite Cluster Evolution",
        fontsize=16,
        fontweight="bold",
        pad=15,
    )
    ax.legend(loc="best", fontsize=12)
    ax.set_xlim(0, STRAIN_MAX * 100)
    plt.tight_layout(pad=2.0)
    pdf.savefig(fig, dpi=300)
    plt.close(fig)

    # GRAPH 7: Nanoscale Lag (Core vs Surface)
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    style_ax(ax)

    ax.plot(
        eps_load * 100,
        smooth(phi_core_load, 7),
        color=C_CORE,
        linewidth=1.8,
        label="Core (Loading)",
        linestyle="-",
    )
    ax.plot(
        eps_load * 100,
        smooth(phi_surf_load, 7),
        color=C_SURF,
        linewidth=1.8,
        label="Surface (Loading)",
        linestyle="--",
    )
    ax.plot(
        eps_unload * 100,
        smooth(phi_core_unload, 7),
        color=C_CORE,
        linewidth=1.8,
        alpha=0.4,
        label="Core (Unloading)",
        linestyle="-",
    )
    ax.plot(
        eps_unload * 100,
        smooth(phi_surf_unload, 7),
        color=C_SURF,
        linewidth=1.8,
        alpha=0.4,
        label="Surface (Unloading)",
        linestyle="--",
    )

    ax.set_xlabel("Compressive Strain, ε (%)", fontsize=14, fontweight="bold")
    ax.set_ylabel("Austenite Fraction, φ$_A$", fontsize=14, fontweight="bold")
    ax.set_title(
        "Nanoscale Lag: Core vs Surface Transformation",
        fontsize=16,
        fontweight="bold",
        pad=15,
    )
    ax.legend(loc="best", fontsize=10, ncol=2)
    ax.set_xlim(0, STRAIN_MAX * 100)
    ax.set_ylim(0, 1.05)
    plt.tight_layout(pad=2.0)
    pdf.savefig(fig, dpi=300)
    plt.close(fig)

    # GRAPH 8: Geometric Breathing (Lattice Distortion)
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    style_ax(ax)

    rg0_load = rg_load[0]
    rg0_unload = rg_unload[0]
    rg_norm_load = (rg_load - rg0_load) / rg0_load * 100
    rg_norm_unload = (rg_unload - rg0_unload) / rg0_unload * 100

    ax.plot(
        eps_load * 100,
        smooth(rg_norm_load, 7),
        color=C_LOAD,
        linewidth=1.8,
        label="Loading",
    )
    ax.plot(
        eps_unload * 100,
        smooth(rg_norm_unload, 7),
        color=C_UNLOAD,
        linewidth=1.8,
        label="Unloading",
    )
    ax.axhline(0, color=C_GRAY, linewidth=1.5, linestyle="--", alpha=0.7)

    ax.set_xlabel("Compressive Strain, ε (%)", fontsize=14, fontweight="bold")
    ax.set_ylabel(
        "Lateral Expansion, ΔR$_g$/R$_g^0$ (%)", fontsize=14, fontweight="bold"
    )
    ax.set_title(
        "Geometric Breathing: Radius of Gyration Change",
        fontsize=16,
        fontweight="bold",
        pad=15,
    )
    ax.legend(loc="best", fontsize=12)
    ax.set_xlim(0, STRAIN_MAX * 100)
    plt.tight_layout(pad=2.0)
    pdf.savefig(fig, dpi=300)
    plt.close(fig)

    # GRAPH 9: Microstructure-Property Coupling Map - SCATTER PLOT VERSION
    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(111)
    style_ax(ax)
    
    # Create smooth curves
    sig_load_smooth = smooth(sig_load, 7)
    sig_unload_smooth = smooth(sig_unload, 7)
    phi_load_smooth = smooth(phi_load, 5)
    phi_unload_smooth = smooth(phi_unload, 5)
    
    # Scatter plot with red-blue colormap
    sc1 = ax.scatter(eps_load*100, sig_load_smooth, c=phi_load_smooth, 
                     cmap='coolwarm', s=35, alpha=0.85, edgecolors='black', 
                     linewidths=0.3, vmin=0.2, vmax=0.8, zorder=3, label='Loading')
    
    sc2 = ax.scatter(eps_unload*100, sig_unload_smooth, c=phi_unload_smooth, 
                     cmap='coolwarm', s=35, alpha=0.6, edgecolors='gray', 
                     linewidths=0.3, vmin=0.2, vmax=0.8, zorder=2, marker='s', label='Unloading')
    
    # Add colorbar
    cbar = fig.colorbar(sc1, ax=ax, pad=0.02, aspect=30)
    cbar.set_label('Austenite Fraction, φ$_A$', fontsize=13, fontweight='bold')
    cbar.ax.tick_params(labelsize=11)
    
    # Create custom legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
               markersize=8, label='Loading', markeredgecolor='black', markeredgewidth=0.5),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='gray', 
               markersize=8, label='Unloading', markeredgecolor='gray', 
               markeredgewidth=0.5, alpha=0.6)
    ]
    ax.legend(handles=legend_elements, loc='best', fontsize=11)
    
    ax.set_xlabel('Compressive Strain, ε (%)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Compressive Stress, |σ$_{zz}$| (GPa)', fontsize=14, fontweight='bold')
    ax.set_title('Microstructure-Property Coupling Map', fontsize=16, fontweight='bold', pad=15)
    ax.set_xlim(0, STRAIN_MAX*100)
    ax.set_ylim(bottom=0)
    plt.tight_layout(pad=2.0)
    pdf.savefig(fig, dpi=300)
    plt.close(fig)

    # Add metadata
    d = pdf.infodict()
    d["Title"] = "NiTi Nanoparticle Superelasticity - Complete Analysis"
    d["Author"] = "MD Simulation"
    d["Subject"] = "Journal-Grade Shape Memory Effect Analysis"
    d["Keywords"] = "NiTi, Superelasticity, Phase Transformation, Nanoparticle"

print(f"\n{'='*80}")
print("ANALYSIS COMPLETE")
print(f"{'='*80}")
print(f"Output Directory: {OUTPUT_DIR}")
print(f"Generated PDF: {pdf_path}")
print(f"Extracted Data: {OUTPUT_DIR}/se_extracted_data.txt")
print(f"{'='*80}\n")

# Print key metrics
sigma_max = float(np.max(np.r_[sig_load, sig_unload]))
phi_min = float(np.min(np.r_[phi_load, phi_unload]))
phi_residual = float(phi_unload[-1])

print("KEY METRICS:")
print(f"  Max stress: {sigma_max:.3f} GPa")
print(f"  Min austenite fraction: {phi_min:.3f}")
print(f"  Residual austenite: {phi_residual:.3f}")
print(f"{'='*80}\n")
