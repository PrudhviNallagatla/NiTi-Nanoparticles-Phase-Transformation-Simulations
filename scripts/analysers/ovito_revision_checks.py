import glob
import json
import math
import warnings
from pathlib import Path

import numpy as np
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier

warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")

BASE = Path(r"d:\MD Sims\SME_SE\data")

CASES = [
    {
        "name": "5nm",
        "raw": BASE / "raw_output" / "niti_5nm_Ni50_single_crystal" / "se",
        "ext": BASE
        / "analysis_output"
        / "niti_5nm_Ni50_single_crystal"
        / "se"
        / "se_extracted_data.txt",
    },
    {
        "name": "10nm",
        "raw": BASE / "raw_output" / "niti_10nm_Ni50_single_crystal" / "se",
        "ext": BASE
        / "analysis_output"
        / "niti_10nm_Ni50_single_crystal"
        / "se"
        / "se_extracted_data.txt",
    },
]

CORE_FRACS = [0.75, 0.80, 0.85]


def load_peak_steps(extracted_file: Path):
    rows = []
    with extracted_file.open("r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            a = line.split()
            eps = float(a[0])
            step = int(float(a[8]))
            rows.append((eps, step))

    peak_steps = [s for eps, s in rows if 0.099 <= eps <= 0.100]
    if not peak_steps:
        peak_steps = [max(rows, key=lambda x: x[0])[1]]
    return sorted(set(peak_steps))


def compute_fractions(data, core_frac: float):
    pos = np.asarray(data.particles.position)
    stype = np.asarray(data.particles["Structure Type"])

    com = pos.mean(axis=0)
    r2 = np.sum((pos - com) ** 2, axis=1)
    rmax = math.sqrt(float(np.max(r2)))

    core = r2 < (core_frac * rmax) ** 2
    surf = ~core

    fcore = float(np.count_nonzero(stype[core] == 3)) / max(int(np.sum(core)), 1)
    fsurf = float(np.count_nonzero(stype[surf] == 3)) / max(int(np.sum(surf)), 1)

    return fcore, fsurf, stype, surf


def main():
    sensitivity_rows = []
    ptm_surface_counts = {}

    for case in CASES:
        peak_steps = load_peak_steps(case["ext"])
        first_dump = glob.glob(str(case["raw"] / "dump.se.*.gz"))[0]
        pipeline = import_file(first_dump)
        pipeline.modifiers.append(
            PolyhedralTemplateMatchingModifier(
                rmsd_cutoff=0.15, output_orientation=False
            )
        )

        for cf in CORE_FRACS:
            fcore_vals = []
            fsurf_vals = []
            for step in peak_steps:
                dump_path = case["raw"] / f"dump.se.{step}.gz"
                pipeline.source.load(str(dump_path))
                data = pipeline.compute()
                fcore, fsurf, stype, surf = compute_fractions(data, cf)
                fcore_vals.append(fcore)
                fsurf_vals.append(fsurf)

                if case["name"] == "10nm" and abs(cf - 0.80) < 1e-9:
                    nonb2 = stype[surf][stype[surf] != 3]
                    vals, counts = np.unique(nonb2, return_counts=True)
                    for v, c in zip(vals, counts):
                        key = str(int(v))
                        ptm_surface_counts[key] = ptm_surface_counts.get(key, 0) + int(
                            c
                        )

            sensitivity_rows.append(
                {
                    "case": case["name"],
                    "core_frac": cf,
                    "fcore_mean": float(np.mean(fcore_vals)),
                    "fsurf_mean": float(np.mean(fsurf_vals)),
                    "delta_mean": float(
                        np.mean(np.array(fcore_vals) - np.array(fsurf_vals))
                    ),
                    "n_frames": len(peak_steps),
                }
            )

    out = {
        "sensitivity": sensitivity_rows,
        "ptm_surface_nonb2_counts_10nm_core080": ptm_surface_counts,
        "ptm_note": "PTM Structure Type 3 corresponds to B2/BCC. Non-3 types group all non-B2 environments and do not uniquely identify B19'.",
    }

    out_path = BASE / "analysis_output" / "uncategorized" / "ovito_revision_checks.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(str(out_path))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
