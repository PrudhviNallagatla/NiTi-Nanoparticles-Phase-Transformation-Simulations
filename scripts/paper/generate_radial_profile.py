#!/usr/bin/env python3
"""
Radial Phase Profile Analysis for NiTi Nanoparticle MD Paper
=============================================================
Computes f_Aus(r/r_max) — the austenite fraction as a function of
normalized radial distance from the nanoparticle center of mass.

This reveals the core-shell structure: high B2 fraction in the interior
decaying toward the disordered surface.

Run from repo root:  python scripts/paper/generate_radial_profile.py

Outputs to: paper/figures/fig_radial_profile.{pdf,png}

Requires: ovito (Python package)
"""

import os
import sys
import numpy as np

try:
    from ovito.io import import_file
    from ovito.modifiers import PolyhedralTemplateMatchingModifier
except ImportError:
    print("ERROR: OVITO Python module not found.")
    print("Install with:  pip install ovito")
    sys.exit(1)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# ==============================================================================
# CONFIGURATION
# ==============================================================================
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUTDIR = os.path.join(REPO_ROOT, "paper", "figures")
os.makedirs(OUTDIR, exist_ok=True)

RAW_BASE = os.path.join(REPO_ROOT, "data", "raw_output")
RMSD = 0.15
N_BINS = 25  # radial bins

# Dumps to analyze
DUMPS = {
    "5nm_eq": {
        "path": os.path.join(
            RAW_BASE, "niti_5nm_Ni50_single_crystal", "se", "dump.se.0.gz"
        ),
        "label": "5 nm, equilibrium",
        "color": "#0072B2",
        "ls": "-",
    },
    "5nm_peak": {
        "path": os.path.join(
            RAW_BASE, "niti_5nm_Ni50_single_crystal", "se", "dump.se.1000000.gz"
        ),
        "label": "5 nm, 10% strain",
        "color": "#0072B2",
        "ls": "--",
    },
    "10nm_eq": {
        "path": os.path.join(
            RAW_BASE, "niti_10nm_Ni50_single_crystal", "se", "dump.se.0.gz"
        ),
        "label": "10 nm, equilibrium",
        "color": "#D55E00",
        "ls": "-",
    },
    "10nm_peak": {
        "path": os.path.join(
            RAW_BASE, "niti_10nm_Ni50_single_crystal", "se", "dump.se.1000000.gz"
        ),
        "label": "10 nm, 10% strain",
        "color": "#D55E00",
        "ls": "--",
    },
}

# ==============================================================================
# STYLE
# ==============================================================================
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif", "Computer Modern Roman"],
        "mathtext.fontset": "cm",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.linewidth": 1.2,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 8.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.minor.size": 3,
        "ytick.minor.size": 3,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        "figure.constrained_layout.use": True,
    }
)

SINGLE_COL = 3.375


# ==============================================================================
# ANALYSIS
# ==============================================================================
def compute_radial_profile(dump_path, n_bins=N_BINS):
    """
    Compute f_Aus as a function of r/r_max.

    Returns:
        bin_centers : (n_bins,) normalized radial positions
        f_aus_bins  : (n_bins,) austenite fraction in each shell
        n_atoms_bins: (n_bins,) atom count in each shell
        n_total     : int, total number of atoms
        r_max       : float, maximum radial distance from COM (Å)
    """
    pipeline = import_file(dump_path)
    ptm = PolyhedralTemplateMatchingModifier(rmsd_cutoff=RMSD, output_orientation=False)
    pipeline.modifiers.append(ptm)
    data = pipeline.compute()

    pos = np.asarray(data.particles.position)
    structs = np.asarray(data.particles["Structure Type"])
    n_total = len(structs)

    # Center of mass
    com = np.mean(pos, axis=0)
    r = np.sqrt(np.sum((pos - com) ** 2, axis=1))
    r_max = np.max(r)

    # Normalize
    r_norm = r / r_max

    # Bin edges (equal spacing in r/r_max from 0 to 1)
    bin_edges = np.linspace(0, 1, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    f_aus_bins = np.zeros(n_bins)
    n_atoms_bins = np.zeros(n_bins)

    for b in range(n_bins):
        mask = (r_norm >= bin_edges[b]) & (r_norm < bin_edges[b + 1])
        n_in_bin = np.sum(mask)
        n_atoms_bins[b] = n_in_bin
        if n_in_bin > 0:
            f_aus_bins[b] = np.count_nonzero(structs[mask] == 3) / n_in_bin

    return bin_centers, f_aus_bins, n_atoms_bins, n_total, r_max


def generate_radial_profile_figure():
    """Generate the radial phase profile figure."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.8), sharey=True)

    # Left panel: equilibrium profiles
    # Right panel: peak strain profiles
    panels = {
        "eq": (ax1, ["5nm_eq", "10nm_eq"], "(a) Equilibrium"),
        "peak": (ax2, ["5nm_peak", "10nm_peak"], r"(b) $\varepsilon = 10\%$"),
    }

    for panel_key, (ax, dump_keys, title) in panels.items():
        for key in dump_keys:
            info = DUMPS[key]
            if not os.path.exists(info["path"]):
                print(f"  WARNING: Missing {info['path']}, skipping {key}")
                continue

            print(f"  Analyzing {key}...")
            r_norm, f_aus, n_atoms, n_total, r_max = compute_radial_profile(
                info["path"]
            )

            ax.plot(
                r_norm,
                f_aus,
                color=info["color"],
                ls=info["ls"],
                lw=1.5,
                label=info["label"],
            )

        # Core fraction boundary line
        ax.axvline(0.8, color="0.6", ls=":", lw=0.8, zorder=0)
        ax.text(
            0.78,
            0.05,
            r"$r_\mathrm{core}$",
            transform=ax.get_xaxis_transform(),
            fontsize=8,
            color="0.5",
            ha="right",
            va="bottom",
        )

        ax.set_xlabel(r"$r / r_\mathrm{max}$")
        ax.set_xlim(0, 1.02)
        ax.set_ylim(-0.02, 1.05)
        ax.set_title(title, fontsize=10, loc="left", fontweight="bold")
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.legend(loc="lower left", frameon=True, framealpha=0.9, edgecolor="0.8")

    ax1.set_ylabel(r"$f_\mathrm{B2}(r)$")

    fig.savefig(os.path.join(OUTDIR, "fig_radial_profile.pdf"))
    fig.savefig(os.path.join(OUTDIR, "fig_radial_profile.png"), dpi=300)
    plt.close(fig)
    print("  [OK] Radial profile figure")


def generate_radial_atom_count():
    """Generate supplementary figure showing atom count per radial bin."""
    fig, ax = plt.subplots(figsize=(SINGLE_COL, 2.5))

    for key in ["5nm_eq", "10nm_eq"]:
        info = DUMPS[key]
        if not os.path.exists(info["path"]):
            continue

        r_norm, _, n_atoms, n_total, r_max = compute_radial_profile(info["path"])

        ax.plot(
            r_norm,
            n_atoms,
            color=info["color"],
            ls="-",
            lw=1.3,
            label=info["label"].split(",")[0],
        )

    ax.set_xlabel(r"$r / r_\mathrm{max}$")
    ax.set_ylabel("Atoms per shell")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(frameon=True, framealpha=0.9, edgecolor="0.8")
    ax.axvline(0.8, color="0.6", ls=":", lw=0.8, zorder=0)

    fig.savefig(os.path.join(OUTDIR, "figS_radial_atom_count.pdf"))
    fig.savefig(os.path.join(OUTDIR, "figS_radial_atom_count.png"), dpi=300)
    plt.close(fig)
    print("  [OK] Radial atom count figure")


# ==============================================================================
# MAIN
# ==============================================================================
if __name__ == "__main__":
    print("=" * 60)
    print("Radial Phase Profile Analysis")
    print(f"Output directory: {OUTDIR}")
    print("=" * 60)

    generate_radial_profile_figure()
    generate_radial_atom_count()

    print("=" * 60)
    print("Done.")
    print("=" * 60)
