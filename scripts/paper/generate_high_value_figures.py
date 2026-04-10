#!/usr/bin/env python3
"""
Generate high-value publication figures for the NiTi/Nitinol nanoparticle paper.

Outputs to: paper/figures/
  - fig8_subcritical_event.{pdf,png}
  - fig9_hysteresis_work.{pdf,png}
  - figS3_se_temperature_control.{pdf,png}
  - figS4_core_surface_deltafb2.{pdf,png}
"""

import os
import json
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUTDIR = os.path.join(REPO_ROOT, "paper", "figures")
os.makedirs(OUTDIR, exist_ok=True)

ANALYSIS = os.path.join(REPO_ROOT, "data", "analysis_output")
RAW = os.path.join(REPO_ROOT, "data", "raw_output")
REVISION_REPORT = os.path.join(
    ANALYSIS,
    "uncategorized",
    "revision_prewrite_analysis.json",
)

PATHS = {
    "5nm": {
        "se": os.path.join(
            ANALYSIS,
            "niti_5nm_Ni50_single_crystal",
            "se",
            "se_extracted_data.txt",
        ),
        "se_raw": os.path.join(
            RAW,
            "niti_5nm_Ni50_single_crystal",
            "se",
            "se_data.txt",
        ),
    },
    "10nm": {
        "se": os.path.join(
            ANALYSIS,
            "niti_10nm_Ni50_single_crystal",
            "se",
            "se_extracted_data.txt",
        ),
        "se_raw": os.path.join(
            RAW,
            "niti_10nm_Ni50_single_crystal",
            "se",
            "se_data.txt",
        ),
        "thermal": os.path.join(
            ANALYSIS,
            "niti_10nm_Ni50_single_crystal",
            "thermal_cycle",
            "extracted_data.txt",
        ),
    },
}

# Protocol landmarks used in existing analysis scripts.
END_STAGE1 = 478571  # end of cooling branch
EQ_OFFSET = 1000
END_LOAD = 1001000

C5 = "#0072B2"
C10 = "#D55E00"
C_CORE = "#009E73"
C_SURF = "#CC79A7"

SINGLE_COL = 3.375
DOUBLE_COL = 7.0

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif", "Computer Modern Roman"],
        "mathtext.fontset": "cm",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 8.5,
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
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        "figure.constrained_layout.use": True,
    }
)


def load_thermal(path):
    d = np.loadtxt(path)
    return {
        "T": d[:, 0],
        "f_aus_core": d[:, 2],
        "clusters": d[:, 7],
        "step": d[:, 9],
    }


def load_se(path):
    d = np.loadtxt(path)
    return {
        "strain": d[:, 0],
        "stress": d[:, 1],
        "f_aus_core": d[:, 3],
        "f_aus_surf": d[:, 4],
        "step": d[:, 8],
    }


def load_se_raw(path):
    # Step T_Global T_Core Lz Strain SigmaZZ(GPa) ...
    d = np.loadtxt(path, comments="#")
    return {
        "step": d[:, 0],
        "T_global": d[:, 1],
        "T_core": d[:, 2],
        "strain": d[:, 4],
    }


def split_load_unload(se_data):
    step = se_data["step"]
    load_mask = (step >= EQ_OFFSET) & (step <= END_LOAD)
    unload_mask = step >= END_LOAD
    load = {k: v[load_mask] for k, v in se_data.items()}
    unload = {k: v[unload_mask] for k, v in se_data.items()}
    return load, unload


def fig8_subcritical_event():
    d = load_thermal(PATHS["10nm"]["thermal"])
    cool = d["step"] <= END_STAGE1

    t = d["T"][cool]
    clusters = d["clusters"][cool]
    fcore = d["f_aus_core"][cool]

    window = (t >= 240.0) & (t <= 300.0)
    t = t[window]
    clusters = clusters[window]
    fcore = fcore[window]

    order = np.argsort(t)
    t = t[order]
    clusters = clusters[order]
    fcore = fcore[order]

    fig, ax1 = plt.subplots(figsize=(SINGLE_COL, 2.8))
    ax2 = ax1.twinx()

    ax1.plot(t, clusters, color=C10, lw=1.6, label="Non-B2 clusters")
    ax2.plot(t, fcore, color=C_CORE, lw=1.4, label=r"$f_{B2,core}$")

    ax1.axvline(272.0, color="0.4", ls="--", lw=0.9)
    ax1.text(272.5, np.max(clusters) * 0.95, "272 K", fontsize=8, color="0.3")

    ax1.set_xlabel("Temperature (K)")
    ax1.set_ylabel("Non-B2 clusters", color=C10)
    ax2.set_ylabel(r"$f_{B2,core}$", color=C_CORE)
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())

    fig.savefig(os.path.join(OUTDIR, "fig8_subcritical_event.pdf"))
    fig.savefig(os.path.join(OUTDIR, "fig8_subcritical_event.png"), dpi=300)
    plt.close(fig)
    print("  [OK] fig8_subcritical_event")


def _hysteresis_area(load, unload):
    eps_l = load["strain"]
    sig_l = load["stress"]
    eps_u = unload["strain"]
    sig_u = unload["stress"]

    order_l = np.argsort(eps_l)
    order_u = np.argsort(eps_u)

    eps_l = eps_l[order_l]
    sig_l = sig_l[order_l]
    eps_u = eps_u[order_u]
    sig_u = sig_u[order_u]

    sig_u_interp = np.interp(eps_l, eps_u, sig_u)
    area = np.trapz(np.abs(sig_l - sig_u_interp), eps_l)
    return eps_l, sig_l, sig_u_interp, area


def fig9_hysteresis_work():
    fig, ax = plt.subplots(figsize=(SINGLE_COL, 2.8))

    reference_areas = {}
    if os.path.exists(REVISION_REPORT):
        with open(REVISION_REPORT, "r", encoding="utf-8") as f:
            report = json.load(f)
        reference_areas = {
            "5nm": report["loop_area_metrics"]["5nm"]["loop_area_gpa"],
            "10nm": report["loop_area_metrics"]["10nm"]["loop_area_gpa"],
        }

    for size, color, label in [("5nm", C5, "5 nm"), ("10nm", C10, "10 nm")]:
        d = load_se(PATHS[size]["se"])
        load, unload = split_load_unload(d)
        eps, sig_l, sig_u_interp, area = _hysteresis_area(load, unload)
        area_to_report = reference_areas.get(size, area)

        ax.plot(eps * 100.0, sig_l, color=color, lw=1.6, ls="-")
        ax.plot(eps * 100.0, sig_u_interp, color=color, lw=1.3, ls="--")
        ax.fill_between(eps * 100.0, sig_l, sig_u_interp, color=color, alpha=0.18)
        ax.text(
            6.0,
            (2.6 if size == "5nm" else 1.35),
            f"{label}: area={area_to_report:.3f} GPa",
            color=color,
            fontsize=8,
        )

    ax.set_xlabel("Engineering strain (%)")
    ax.set_ylabel("Effective virial stress (GPa)")
    ax.set_xlim(0.0, 10.0)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    fig.savefig(os.path.join(OUTDIR, "fig9_hysteresis_work.pdf"))
    fig.savefig(os.path.join(OUTDIR, "fig9_hysteresis_work.png"), dpi=300)
    plt.close(fig)
    print("  [OK] fig9_hysteresis_work")


def figS3_se_temperature_control():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 2.8), sharey=True)

    for ax, size, color in [(ax1, "5nm", C5), (ax2, "10nm", C10)]:
        d = load_se_raw(PATHS[size]["se_raw"])
        step = d["step"]
        mask = step >= EQ_OFFSET

        eps = d["strain"][mask] * 100.0
        tg = d["T_global"][mask]
        tc = d["T_core"][mask]

        ax.plot(eps, tg, color=color, lw=1.4, ls="-", label="Global")
        ax.plot(eps, tc, color=color, lw=1.2, ls="--", label="Core")
        ax.axhline(500.0, color="0.5", ls=":", lw=0.8)
        ax.set_title(f"{size.replace('nm', ' nm')}")
        ax.set_xlabel("Engineering strain (%)")
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.legend(frameon=True, framealpha=0.9, edgecolor="0.8")

    ax1.set_ylabel("Temperature (K)")

    fig.savefig(os.path.join(OUTDIR, "figS3_se_temperature_control.pdf"))
    fig.savefig(os.path.join(OUTDIR, "figS3_se_temperature_control.png"), dpi=300)
    plt.close(fig)
    print("  [OK] figS3_se_temperature_control")


def figS4_core_surface_deltafb2():
    fig, ax = plt.subplots(figsize=(SINGLE_COL, 2.8))

    for size, color, label in [("5nm", C5, "5 nm"), ("10nm", C10, "10 nm")]:
        d = load_se(PATHS[size]["se"])
        load, unload = split_load_unload(d)

        step = np.concatenate([load["step"], unload["step"]])
        strain = np.concatenate([load["strain"], unload["strain"]]) * 100.0
        delta = np.concatenate(
            [
                load["f_aus_core"] - load["f_aus_surf"],
                unload["f_aus_core"] - unload["f_aus_surf"],
            ]
        )

        order = np.argsort(step)
        ax.plot(strain[order], delta[order], color=color, lw=1.5, label=label)

    ax.axhline(0.0, color="0.5", ls=":", lw=0.8)
    ax.set_xlabel("Engineering strain (%)")
    ax.set_ylabel(r"$\Delta f_{B2} = f_{B2,core} - f_{B2,surf}$")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(frameon=True, framealpha=0.9, edgecolor="0.8")

    fig.savefig(os.path.join(OUTDIR, "figS4_core_surface_deltafb2.pdf"))
    fig.savefig(os.path.join(OUTDIR, "figS4_core_surface_deltafb2.png"), dpi=300)
    plt.close(fig)
    print("  [OK] figS4_core_surface_deltafb2")


if __name__ == "__main__":
    print("=" * 60)
    print("Generating high-value figures")
    print(f"Output directory: {OUTDIR}")
    print("=" * 60)

    fig8_subcritical_event()
    fig9_hysteresis_work()
    figS3_se_temperature_control()
    figS4_core_surface_deltafb2()

    print("=" * 60)
    print("Done")
    print("=" * 60)
