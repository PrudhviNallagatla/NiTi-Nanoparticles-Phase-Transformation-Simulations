#!/usr/bin/env python3
"""
Publication-Quality Figure Generator for NiTi Nanoparticle MD Paper
===================================================================
Generates all matplotlib-based figures from OVITO-extracted data.
Run from repo root:  python scripts/paper/generate_data_figures.py

Outputs to: paper/figures/

Figure Plan:
  Fig 2  - Thermal cycle: f_Aus(T) hysteresis for 5 nm and 10 nm
  Fig 3  - Thermal cycle: PE_core(T) and Rg(T) dual panel
  Fig 4  - SE: stress-strain curves (5 nm and 10 nm overlay)
  Fig 5  - SE: f_Aus(strain) decomposed into core/surface/total
  Fig S1 - SE: cluster count vs strain
  Fig S2 - Thermal: core vs surface austenite fraction

"""

import os
import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

# PATH CONFIGURATION
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

DATA = {
    "5nm": {
        "thermal": os.path.join(
            REPO_ROOT,
            "data",
            "analysis_output",
            "niti_5nm_Ni50_single_crystal",
            "thermal_cycle",
            "extracted_data.txt",
        ),
        "se": os.path.join(
            REPO_ROOT,
            "data",
            "analysis_output",
            "niti_5nm_Ni50_single_crystal",
            "se",
            "se_extracted_data.txt",
        ),
        "label": "d = 5 nm",
        "N": 4745,
        "N_core": 2513,
    },
    "10nm": {
        "thermal": os.path.join(
            REPO_ROOT,
            "data",
            "analysis_output",
            "niti_10nm_Ni50_single_crystal",
            "thermal_cycle",
            "extracted_data.txt",
        ),
        "se": os.path.join(
            REPO_ROOT,
            "data",
            "analysis_output",
            "niti_10nm_Ni50_single_crystal",
            "se",
            "se_extracted_data.txt",
        ),
        "label": "d = 10 nm",
        "N": 37869,
        "N_core": 20261,
    },
}

OUTDIR = os.path.join(REPO_ROOT, "paper", "figures")
os.makedirs(OUTDIR, exist_ok=True)

# STYLE CONFIGURATION - publication quality
# Colors
C5 = "#0072B2"  # blue  (5 nm)
C10 = "#D55E00"  # vermillion (10 nm)
C5_LIGHT = "#56B4E9"
C10_LIGHT = "#E69F00"
C_CORE = "#009E73"  # green
C_SURF = "#CC79A7"  # pink
C_TOTAL = "#444444"  # dark gray

# Matplotlib rc
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif", "Computer Modern Roman"],
        "mathtext.fontset": "cm",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.minor.width": 0.6,
        "ytick.minor.width": 0.6,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.minor.size": 3,
        "ytick.minor.size": 3,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "lines.linewidth": 1.5,
        "lines.markersize": 3,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        "figure.constrained_layout.use": True,
    }
)

SINGLE_COL = 3.375  # inches (APS/Elsevier single column)
DOUBLE_COL = 7.0  # inches (double column)


# DATA LOADING
def load_thermal(path):
    """Load thermal cycle extracted data.
    Columns: Temp f_Aus f_Aus_Core f_Aus_Surf PE_Core Shear Stress Clusters Rg Step
    """
    d = np.loadtxt(path)
    return {
        "T": d[:, 0],
        "f_aus": d[:, 1],
        "f_aus_core": d[:, 2],
        "f_aus_surf": d[:, 3],
        "pe_core": d[:, 4],
        "shear": d[:, 5],
        "stress": d[:, 6],
        "clusters": d[:, 7],
        "rg": d[:, 8],
        "step": d[:, 9],
    }


def load_se(path):
    """Load SE extracted data.
    Columns: Strain Stress_GPa f_Aus f_Aus_Core f_Aus_Surf PE_Core Clusters Rg Step
    """
    d = np.loadtxt(path)
    return {
        "strain": d[:, 0],
        "stress": d[:, 1],
        "f_aus": d[:, 2],
        "f_aus_core": d[:, 3],
        "f_aus_surf": d[:, 4],
        "pe_core": d[:, 5],
        "clusters": d[:, 6],
        "rg": d[:, 7],
        "step": d[:, 8],
    }


def smooth(y, sigma=3):
    """Gaussian smoothing for noisy MD data (NumPy-only implementation)."""
    y = np.asarray(y, dtype=float)
    if sigma <= 0:
        return y
    radius = int(max(1, round(3 * sigma)))
    x = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-(x**2) / (2.0 * sigma**2))
    kernel /= np.sum(kernel)
    return np.convolve(y, kernel, mode="same")


def split_cool_heat(data):
    """Split thermal cycle into cooling (500->200 K) and heating (200->500 K) branches."""
    step = data["step"]
    # Protocol: hold 500K -> cool to 200K -> hold 200K -> heat to 500K
    # From thermal_cycle.lmp: RUN_HOLD=50000, RUN_RAMP = (500-200)/0.35/0.002 = 428571
    # END_STAGE0 = 50000, END_STAGE1 = 478571, END_STAGE2 = 528571, END_STAGE3 = 957142
    END_STAGE0 = 50000
    END_STAGE1 = 478571
    END_STAGE2 = 528571
    END_STAGE3 = 957142

    # Cooling: stage 0 (hold at 500K) + stage 1 (cool 500->200K)
    cool_mask = step <= END_STAGE1
    # Heating: stage 2 (hold at 200K) + stage 3 (heat 200->500K)
    heat_mask = step >= END_STAGE2

    cool = {k: v[cool_mask] for k, v in data.items()}
    heat = {k: v[heat_mask] for k, v in data.items()}
    return cool, heat


def split_load_unload(data):
    """Split SE data into loading and unloading branches."""
    step = data["step"]
    EQ_OFFSET = 1000
    END_LOAD = 1001000
    load_mask = (step >= EQ_OFFSET) & (step <= END_LOAD)
    unload_mask = step >= END_LOAD
    load = {k: v[load_mask] for k, v in data.items()}
    unload = {k: v[unload_mask] for k, v in data.items()}
    return load, unload


# FIGURE 2: f_Aus(T) — Thermal Hysteresis
def fig2_thermal_hysteresis():
    """f_Aus vs Temperature for both NP sizes."""
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(SINGLE_COL, 4.8), sharex=True, sharey=True
    )

    for ax, size, color, c_light in [
        (ax1, "5nm", C5, C5_LIGHT),
        (ax2, "10nm", C10, C10_LIGHT),
    ]:
        d = load_thermal(DATA[size]["thermal"])
        cool, heat = split_cool_heat(d)

        # Raw data as faint scatter
        ax.scatter(
            cool["T"],
            cool["f_aus"],
            s=1,
            alpha=0.15,
            color=c_light,
            rasterized=True,
            zorder=1,
        )
        ax.scatter(
            heat["T"],
            heat["f_aus"],
            s=1,
            alpha=0.15,
            color=c_light,
            rasterized=True,
            zorder=1,
        )

        # Smoothed trend lines
        ax.plot(
            cool["T"],
            smooth(cool["f_aus"], sigma=5),
            color=color,
            ls="-",
            lw=1.8,
            label="Cooling",
            zorder=3,
        )
        ax.plot(
            heat["T"],
            smooth(heat["f_aus"], sigma=5),
            color=color,
            ls="--",
            lw=1.8,
            label="Heating",
            zorder=3,
        )

        ax.set_xlabel(r"Temperature (K)")
        ax.set_title(DATA[size]["label"], fontsize=10, fontweight="bold")
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_xlim(190, 510)
        ax.legend(loc="lower left", frameon=True, framealpha=0.9, edgecolor="0.8")

    ax1.set_ylabel(r"$f_\mathrm{B2}$")
    ax1.text(0.03, 0.95, "(a)", transform=ax1.transAxes, fontweight="bold", va="top")
    ax2.text(0.03, 0.95, "(b)", transform=ax2.transAxes, fontweight="bold", va="top")

    fig.savefig(os.path.join(OUTDIR, "fig2_thermal_hysteresis.pdf"))
    fig.savefig(os.path.join(OUTDIR, "fig2_thermal_hysteresis.png"), dpi=300)
    plt.close(fig)
    print("  [OK] Fig 2: Thermal hysteresis")


# FIGURE 3: PE_core(T) and Rg(T)
def fig3_pe_rg_thermal():
    """Dual-panel: PE_core and Rg vs Temperature."""
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COL, 4.5), sharex="col")

    for j, (size, color, c_light) in enumerate(
        [("5nm", C5, C5_LIGHT), ("10nm", C10, C10_LIGHT)]
    ):
        d = load_thermal(DATA[size]["thermal"])
        cool, heat = split_cool_heat(d)

        # --- PE_core ---
        ax = axes[0, j]
        ax.scatter(
            cool["T"], cool["pe_core"], s=1, alpha=0.12, color=c_light, rasterized=True
        )
        ax.scatter(
            heat["T"], heat["pe_core"], s=1, alpha=0.12, color=c_light, rasterized=True
        )
        ax.plot(
            cool["T"], smooth(cool["pe_core"], 5), color=color, ls="-", label="Cooling"
        )
        ax.plot(
            heat["T"], smooth(heat["pe_core"], 5), color=color, ls="--", label="Heating"
        )
        if j == 0:
            ax.set_ylabel(r"$\langle E_\mathrm{pot}\rangle_\mathrm{core}$ (eV/atom)")
        ax.set_title(DATA[size]["label"], fontsize=10, fontweight="bold")
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.legend(
            loc="upper left", frameon=True, framealpha=0.9, edgecolor="0.8", fontsize=8
        )

        # --- Rg ---
        ax = axes[1, j]
        ax.scatter(
            cool["T"], cool["rg"], s=1, alpha=0.12, color=c_light, rasterized=True
        )
        ax.scatter(
            heat["T"], heat["rg"], s=1, alpha=0.12, color=c_light, rasterized=True
        )
        ax.plot(cool["T"], smooth(cool["rg"], 5), color=color, ls="-", label="Cooling")
        ax.plot(heat["T"], smooth(heat["rg"], 5), color=color, ls="--", label="Heating")
        if j == 0:
            ax.set_ylabel(r"$R_\mathrm{g}$ (Å)")
        ax.set_xlabel(r"Temperature (K)")
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

    labels = ["(a)", "(b)", "(c)", "(d)"]
    for i, ax in enumerate(axes.flat):
        ax.text(
            0.03, 0.95, labels[i], transform=ax.transAxes, fontweight="bold", va="top"
        )

    fig.savefig(os.path.join(OUTDIR, "fig3_pe_rg_thermal.pdf"))
    fig.savefig(os.path.join(OUTDIR, "fig3_pe_rg_thermal.png"), dpi=300)
    plt.close(fig)
    print("  [OK] Fig 3: PE and Rg thermal")


# FIGURE 4: Stress-Strain (Superelasticity)
def fig4_stress_strain():
    """Compressive stress-strain curves for both NP sizes."""
    fig, ax = plt.subplots(figsize=(SINGLE_COL, 2.8))

    for size, color, marker in [("5nm", C5, ""), ("10nm", C10, "")]:
        d = load_se(DATA[size]["se"])
        load, unload = split_load_unload(d)

        # Smooth the noisy virial stress
        strain_load = load["strain"] * 100  # percent
        stress_load = smooth(load["stress"], sigma=4)
        strain_unload = unload["strain"] * 100
        stress_unload = smooth(unload["stress"], sigma=4)

        ax.plot(
            strain_load,
            stress_load,
            color=color,
            ls="-",
            lw=1.5,
            label=DATA[size]["label"] + " load",
            zorder=3,
        )
        ax.plot(
            strain_unload,
            stress_unload,
            color=color,
            ls="--",
            lw=1.3,
            label=DATA[size]["label"] + " unload",
            zorder=3,
        )

    ax.set_xlabel("Engineering strain (%)")
    ax.set_ylabel(r"Effective virial stress (GPa)")
    ax.set_xlim(-0.2, 10.5)
    ax.axhline(0, color="0.7", lw=0.5, zorder=0)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(
        loc="upper left", frameon=True, framealpha=0.9, edgecolor="0.8", fontsize=8
    )

    fig.savefig(os.path.join(OUTDIR, "fig4_stress_strain.pdf"))
    fig.savefig(os.path.join(OUTDIR, "fig4_stress_strain.png"), dpi=300)
    plt.close(fig)
    print("  [OK] Fig 4: Stress-strain")


# FIGURE 5: f_Aus decomposition during SE
def fig5_phase_fraction_se():
    """f_Aus (total, core, surface) vs strain for both sizes."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 2.8), sharey=True)

    for ax, size in [(ax1, "5nm"), (ax2, "10nm")]:
        d = load_se(DATA[size]["se"])
        load, unload = split_load_unload(d)

        strain_all = np.concatenate([load["strain"], unload["strain"]]) * 100
        f_total = np.concatenate([load["f_aus"], unload["f_aus"]])
        f_core = np.concatenate([load["f_aus_core"], unload["f_aus_core"]])
        f_surf = np.concatenate([load["f_aus_surf"], unload["f_aus_surf"]])

        # Sort by step to get continuous path
        steps = np.concatenate([load["step"], unload["step"]])
        order = np.argsort(steps)
        strain_all = strain_all[order]
        f_total = f_total[order]
        f_core = f_core[order]
        f_surf = f_surf[order]

        # Plot with smoothing
        ax.plot(
            strain_all, smooth(f_total, 4), color=C_TOTAL, ls="-", lw=1.5, label="Total"
        )
        ax.plot(
            strain_all, smooth(f_core, 4), color=C_CORE, ls="-", lw=1.3, label="Core"
        )
        ax.plot(
            strain_all, smooth(f_surf, 4), color=C_SURF, ls="-", lw=1.3, label="Surface"
        )

        # Add arrows to show load/unload direction
        n = len(strain_all)
        mid_load = n // 5
        mid_unload = 3 * n // 4
        for idx in [mid_load, mid_unload]:
            if idx + 1 < n:
                dx = strain_all[idx + 1] - strain_all[idx]
                dy = smooth(f_total, 4)[idx + 1] - smooth(f_total, 4)[idx]
                ax.annotate(
                    "",
                    xy=(strain_all[idx] + dx, smooth(f_total, 4)[idx] + dy),
                    xytext=(strain_all[idx], smooth(f_total, 4)[idx]),
                    arrowprops=dict(arrowstyle="->", color=C_TOTAL, lw=1.5),
                )

        ax.set_xlabel("Engineering strain (%)")
        ax.set_title(DATA[size]["label"], fontsize=10, fontweight="bold")
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.legend(
            loc="lower left", frameon=True, framealpha=0.9, edgecolor="0.8", fontsize=8
        )

    ax1.set_ylabel(r"$f_\mathrm{B2}$")
    ax1.text(0.03, 0.95, "(a)", transform=ax1.transAxes, fontweight="bold", va="top")
    ax2.text(0.03, 0.95, "(b)", transform=ax2.transAxes, fontweight="bold", va="top")

    fig.savefig(os.path.join(OUTDIR, "fig5_phase_fraction_se.pdf"))
    fig.savefig(os.path.join(OUTDIR, "fig5_phase_fraction_se.png"), dpi=300)
    plt.close(fig)
    print("  [OK] Fig 5: Phase fraction SE")


# FIGURE S1: Cluster count evolution during SE
def figS1_clusters_se():
    """Number of non-B2 clusters vs strain."""
    fig, ax = plt.subplots(figsize=(SINGLE_COL, 2.8))

    for size, color in [("5nm", C5), ("10nm", C10)]:
        d = load_se(DATA[size]["se"])
        load, unload = split_load_unload(d)

        strain_all = np.concatenate([load["strain"], unload["strain"]]) * 100
        clust = np.concatenate([load["clusters"], unload["clusters"]])
        steps = np.concatenate([load["step"], unload["step"]])
        order = np.argsort(steps)

        ax.plot(
            strain_all[order],
            smooth(clust[order], 3),
            color=color,
            lw=1.3,
            label=DATA[size]["label"],
        )

    ax.set_xlabel("Engineering strain (%)")
    ax.set_ylabel("Non-B2 cluster count")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(frameon=True, framealpha=0.9, edgecolor="0.8")

    fig.savefig(os.path.join(OUTDIR, "figS1_clusters_se.pdf"))
    fig.savefig(os.path.join(OUTDIR, "figS1_clusters_se.png"), dpi=300)
    plt.close(fig)
    print("  [OK] Fig S1: Cluster count SE")


# FIGURE S2: Core vs Surface in thermal cycle
def figS2_core_surf_thermal():
    """f_Aus_core and f_Aus_surf vs T for both sizes."""
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(SINGLE_COL, 4.8), sharex=True, sharey=True
    )

    for ax, size in [(ax1, "5nm"), (ax2, "10nm")]:
        d = load_thermal(DATA[size]["thermal"])
        cool, heat = split_cool_heat(d)

        for branch, ls in [(cool, "-"), (heat, "--")]:
            ax.plot(
                branch["T"],
                smooth(branch["f_aus_core"], 5),
                color=C_CORE,
                ls=ls,
                lw=1.3,
                label="Core" if ls == "-" else None,
            )
            ax.plot(
                branch["T"],
                smooth(branch["f_aus_surf"], 5),
                color=C_SURF,
                ls=ls,
                lw=1.3,
                label="Surface" if ls == "-" else None,
            )

        ax.set_xlabel("Temperature (K)")
        ax.set_title(DATA[size]["label"], fontsize=10, fontweight="bold")
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.legend(frameon=True, framealpha=0.9, edgecolor="0.8", fontsize=8)

    ax1.set_ylabel(r"$f_\mathrm{B2}$")
    ax1.text(0.03, 0.95, "(a)", transform=ax1.transAxes, fontweight="bold", va="top")
    ax2.text(0.03, 0.95, "(b)", transform=ax2.transAxes, fontweight="bold", va="top")

    fig.savefig(os.path.join(OUTDIR, "figS2_core_surf_thermal.pdf"))
    fig.savefig(os.path.join(OUTDIR, "figS2_core_surf_thermal.png"), dpi=300)
    plt.close(fig)
    print("  [OK] Fig S2: Core vs Surface thermal")


# FIGURE 6 (bonus): Size-comparison summary panel
def fig6_size_comparison():
    """4-panel: stress-strain overlay, f_aus overlay, PE_core(strain), Rg(strain)."""
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COL, 4.5))

    for size, color, c_light in [("5nm", C5, C5_LIGHT), ("10nm", C10, C10_LIGHT)]:
        d = load_se(DATA[size]["se"])
        load, unload = split_load_unload(d)
        strain_all = np.concatenate([load["strain"], unload["strain"]]) * 100
        steps = np.concatenate([load["step"], unload["step"]])
        order = np.argsort(steps)
        strain_sorted = strain_all[order]

        # (a) stress-strain
        ax = axes[0, 0]
        ax.plot(
            load["strain"] * 100,
            smooth(load["stress"], 4),
            color=color,
            ls="-",
            label=DATA[size]["label"],
        )
        ax.plot(
            unload["strain"] * 100, smooth(unload["stress"], 4), color=color, ls="--"
        )

        # (b) f_aus(strain)
        ax = axes[0, 1]
        f_tot = np.concatenate([load["f_aus"], unload["f_aus"]])
        ax.plot(
            strain_sorted,
            smooth(f_tot[order], 4),
            color=color,
            label=DATA[size]["label"],
        )

        # (c) PE_core(strain)
        ax = axes[1, 0]
        pe = np.concatenate([load["pe_core"], unload["pe_core"]])
        ax.plot(
            strain_sorted, smooth(pe[order], 4), color=color, label=DATA[size]["label"]
        )

        # (d) Rg(strain)
        ax = axes[1, 1]
        rg = np.concatenate([load["rg"], unload["rg"]])
        ax.plot(
            strain_sorted, smooth(rg[order], 4), color=color, label=DATA[size]["label"]
        )

    ylabels = [
        r"Effective virial stress (GPa)",
        r"$f_\mathrm{B2}$",
        r"$\langle E_\mathrm{pot}\rangle_\mathrm{core}$ (eV/atom)",
        r"$R_\mathrm{g}$ (Å)",
    ]
    for i, ax in enumerate(axes.flat):
        ax.set_xlabel("Engineering strain (%)")
        ax.set_ylabel(ylabels[i])
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.legend(frameon=True, framealpha=0.9, edgecolor="0.8", fontsize=8)
        ax.text(
            0.03,
            0.95,
            f"({'abcd'[i]})",
            transform=ax.transAxes,
            fontweight="bold",
            va="top",
        )

    fig.savefig(os.path.join(OUTDIR, "fig6_size_comparison.pdf"))
    fig.savefig(os.path.join(OUTDIR, "fig6_size_comparison.png"), dpi=300)
    plt.close(fig)
    print("  [OK] Fig 6: Size comparison")


# MAIN
if __name__ == "__main__":
    print("=" * 60)
    print("Generating publication figures...")
    print(f"Output directory: {OUTDIR}")
    print("=" * 60)

    # Verify data exists
    for size, info in DATA.items():
        for key in ["thermal", "se"]:
            if not os.path.exists(info[key]):
                print(f"ERROR: Missing {key} data for {size}: {info[key]}")
                sys.exit(1)

    fig2_thermal_hysteresis()
    fig3_pe_rg_thermal()
    fig4_stress_strain()
    fig5_phase_fraction_se()
    fig6_size_comparison()
    figS1_clusters_se()
    figS2_core_surf_thermal()

    print("=" * 60)
    print(f"All figures saved to: {OUTDIR}")
    print("=" * 60)
