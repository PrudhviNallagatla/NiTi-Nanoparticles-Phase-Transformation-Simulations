import json
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parents[2]
DATA = BASE / "data" / "analysis_output" / "uncategorized"
OUT = DATA / "rphase_proxy_metrics.json"


CASES = {
    "5nm": {
        "rdf": DATA / "rdf_nonb2_peak_5nm.csv",
        "angles": DATA / "bond_angle_nonb2_peak_5nm.csv",
    },
    "10nm": {
        "rdf": DATA / "rdf_nonb2_peak_10nm.csv",
        "angles": DATA / "bond_angle_nonb2_peak_10nm.csv",
    },
}


def integrate_window(x: np.ndarray, y: np.ndarray, lo: float, hi: float) -> float:
    mask = (x >= lo) & (x <= hi)
    if not np.any(mask):
        return 0.0
    return float(np.trapz(y[mask], x[mask]))


def analyze_case(rdf_csv: Path, angle_csv: Path) -> dict:
    rdf = np.loadtxt(rdf_csv, delimiter=",", skiprows=1)
    ang = np.loadtxt(angle_csv, delimiter=",", skiprows=1)

    r = rdf[:, 0]
    g = rdf[:, 1]
    theta = ang[:, 0]
    pdf = ang[:, 1]

    peak_idx = int(np.argmax(g))
    first_peak_r = float(r[peak_idx])

    # RDF proxies discussed in manuscript/review comments.
    rdf_26_28 = integrate_window(r, g, 2.6, 2.8)
    rdf_28_30 = integrate_window(r, g, 2.8, 3.0)
    rdf_40_45 = integrate_window(r, g, 4.0, 4.5)
    rdf_45_50 = integrate_window(r, g, 4.5, 5.0)

    first_shell_shoulder_ratio = rdf_26_28 / max(rdf_28_30, 1e-12)
    second_neighbor_symmetry_ratio = rdf_40_45 / max(rdf_45_50, 1e-12)

    # Bond-angle proxies: 58-62 (shared), 90-95 (B19'-leaning), 118-122 (R-leaning).
    ang_58_62 = integrate_window(theta, pdf, 58.0, 62.0)
    ang_90_95 = integrate_window(theta, pdf, 90.0, 95.0)
    ang_118_122 = integrate_window(theta, pdf, 118.0, 122.0)

    ratio_90_to_60 = ang_90_95 / max(ang_58_62, 1e-12)
    ratio_120_to_60 = ang_118_122 / max(ang_58_62, 1e-12)

    top_idx = np.argsort(pdf)[-8:][::-1]
    top_angles = [{"theta_deg": float(theta[i]), "pdf": float(pdf[i])} for i in top_idx]

    if ratio_120_to_60 > ratio_90_to_60:
        verdict = (
            "R-phase-leaning proxy signature (120 deg contribution exceeds 90-95 deg), "
            "but non-unique without explicit template classification"
        )
    else:
        verdict = "B19'-leaning proxy signature, but non-unique without explicit template classification"

    return {
        "first_peak_r_angstrom": first_peak_r,
        "rdf_window_integrals": {
            "2.6_2.8": rdf_26_28,
            "2.8_3.0": rdf_28_30,
            "4.0_4.5": rdf_40_45,
            "4.5_5.0": rdf_45_50,
        },
        "rdf_proxy_ratios": {
            "first_shell_shoulder_ratio_2p6_2p8_over_2p8_3p0": first_shell_shoulder_ratio,
            "second_neighbor_symmetry_ratio_4p0_4p5_over_4p5_5p0": second_neighbor_symmetry_ratio,
        },
        "angle_window_integrals": {
            "58_62": ang_58_62,
            "90_95": ang_90_95,
            "118_122": ang_118_122,
        },
        "angle_proxy_ratios": {
            "ratio_90_95_to_58_62": ratio_90_to_60,
            "ratio_118_122_to_58_62": ratio_120_to_60,
        },
        "top_angle_peaks": top_angles,
        "proxy_verdict": verdict,
    }


def main() -> None:
    results = {}
    for name, paths in CASES.items():
        if not paths["rdf"].exists() or not paths["angles"].exists():
            raise FileNotFoundError(
                f"Missing inputs for {name}: {paths['rdf']} or {paths['angles']}"
            )
        results[name] = analyze_case(paths["rdf"], paths["angles"])

    out = {
        "method": (
            "Proxy discrimination from existing non-B2 RDF and bond-angle distributions. "
            "This does not replace explicit PTM/CNA with an R-phase template."
        ),
        "cases": results,
    }

    OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(str(OUT))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
