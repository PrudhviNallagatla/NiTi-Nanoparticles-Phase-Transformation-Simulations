from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "analysis_output" / "uncategorized"
FIG = ROOT / "paper" / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def plot_rmsd(report: dict):
    rows = []
    for size in ["5nm", "10nm"]:
        core = report["ptm_rmsd_core_surface_equilibrium"][size]["core"]
        surf = report["ptm_rmsd_core_surface_equilibrium"][size]["surface"]
        rows.append((size, "Core", core["median"], core["p90"], core["p99"]))
        rows.append((size, "Surface", surf["median"], surf["p90"], surf["p99"]))

    x = np.arange(len(rows))
    med = [r[2] for r in rows]
    p90 = [r[3] for r in rows]
    p99 = [r[4] for r in rows]
    labels = [f"{r[0]}\n{r[1]}" for r in rows]

    plt.figure(figsize=(8, 4.5))
    plt.plot(x, med, marker="o", label="Median")
    plt.plot(x, p90, marker="s", label="P90")
    plt.plot(x, p99, marker="^", label="P99")
    plt.xticks(x, labels)
    plt.ylabel("PTM RMSD")
    plt.title("Core vs surface PTM RMSD summary at equilibrium")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "figS5_ptm_rmsd_core_surface.pdf")
    plt.close()


def plot_rdf():
    d5 = pd.read_csv(DATA / "rdf_nonb2_peak_5nm.csv")
    d10 = pd.read_csv(DATA / "rdf_nonb2_peak_10nm.csv")

    plt.figure(figsize=(7, 4.5))
    plt.plot(d5["r_angstrom"], d5["g_r"], label="5 nm non-B2")
    plt.plot(d10["r_angstrom"], d10["g_r"], label="10 nm non-B2")
    plt.xlabel("r (A)")
    plt.ylabel("g(r)")
    plt.title("Peak-strain non-B2 RDF")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "figS6_rdf_nonb2_peak.pdf")
    plt.close()


def plot_bond_angle():
    d5 = pd.read_csv(DATA / "bond_angle_nonb2_peak_5nm.csv")
    d10 = pd.read_csv(DATA / "bond_angle_nonb2_peak_10nm.csv")

    plt.figure(figsize=(7, 4.5))
    plt.plot(d5["theta_deg"], d5["pdf"], label="5 nm non-B2")
    plt.plot(d10["theta_deg"], d10["pdf"], label="10 nm non-B2")
    plt.xlabel("Bond angle (deg)")
    plt.ylabel("PDF")
    plt.title("Peak-strain non-B2 bond-angle distribution")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "figS7_bond_angle_nonb2_peak.pdf")
    plt.close()


def plot_nucleation_tracking():
    d5 = pd.read_csv(DATA / "nucleation_tracking_5nm.csv")
    d10 = pd.read_csv(DATA / "nucleation_tracking_10nm.csv")

    fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    ax[0].plot(d5["strain"], d5["cluster_centroid_r_over_rmax"], lw=1.5)
    ax[0].set_title("5 nm")
    ax[0].set_xlabel("Strain")
    ax[0].set_ylabel("Cluster centroid r/rmax")
    ax[0].grid(alpha=0.25)

    ax[1].plot(d10["strain"], d10["cluster_centroid_r_over_rmax"], lw=1.5)
    ax[1].set_title("10 nm")
    ax[1].set_xlabel("Strain")
    ax[1].grid(alpha=0.25)

    fig.suptitle("Largest core non-B2 cluster radial position during loading")
    plt.tight_layout()
    plt.savefig(FIG / "figS8_nucleation_tracking.pdf")
    plt.close()


def plot_rphase_panel(report: dict):
    templates = report["ptm_r_phase_support"]["available_templates"]
    keyframes = report["ptm_keyframe_phase_fractions"]

    labels = []
    bcc = []
    other = []
    for size in ["5nm", "10nm"]:
        for frame in [
            "equilibrium_500K",
            "cooling_near_272K",
            "peak_strain",
            "recovered",
        ]:
            entry = keyframes[size][frame]["fractions"]
            labels.append(f"{size}\n{frame.split('_')[0]}")
            bcc.append(entry["BCC"]["fraction"])
            other.append(entry["Other"]["fraction"])

    x = np.arange(len(labels))
    width = 0.38

    plt.figure(figsize=(11, 4.5))
    plt.bar(x - width / 2, bcc, width=width, label="PTM BCC")
    plt.bar(x + width / 2, other, width=width, label="PTM Other")
    plt.xticks(x, labels, rotation=30, ha="right")
    plt.ylabel("Fraction")
    plt.title("PTM key-frame template fractions (no explicit R-phase template)")
    plt.text(
        0.01,
        0.02,
        f"Templates: {', '.join(templates)}",
        transform=plt.gca().transAxes,
        fontsize=8,
    )
    plt.grid(axis="y", alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "figS9_ptm_keyframe_templates.pdf")
    plt.close()


def plot_loop_panel(report: dict):
    m5 = report["loop_area_metrics"]["5nm"]
    m10 = report["loop_area_metrics"]["10nm"]

    labels = ["5 nm", "10 nm"]
    raw = [m5["loop_area_gpa"], m10["loop_area_gpa"]]
    norm = [m5["normalized_loop_area"], m10["normalized_loop_area"]]

    x = np.arange(len(labels))
    width = 0.35

    plt.figure(figsize=(7, 4.5))
    plt.bar(x - width / 2, raw, width=width, label="Raw loop area (GPa)")
    plt.bar(x + width / 2, norm, width=width, label="Normalized loop area")
    plt.xticks(x, labels)
    plt.title("Raw and normalized hysteresis loop metrics")
    plt.grid(axis="y", alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "figS10_loop_area_normalized.pdf")
    plt.close()


def main():
    report = load_json(DATA / "revision_prewrite_analysis.json")
    plot_rmsd(report)
    plot_rdf()
    plot_bond_angle()
    plot_nucleation_tracking()
    plot_rphase_panel(report)
    plot_loop_panel(report)
    print("Generated revision supplementary figures in", FIG)


if __name__ == "__main__":
    main()
