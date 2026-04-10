import json
from pathlib import Path

import numpy as np
from ovito.io import import_file
from ovito.modifiers import (
    AcklandJonesModifier,
    CommonNeighborAnalysisModifier,
    PolyhedralTemplateMatchingModifier,
)

BASE = Path(__file__).resolve().parents[2]
DATA = BASE / "data"
OUT = DATA / "analysis_output" / "uncategorized" / "rphase_cna_aj_scan.json"

CASES = {
    "5nm": {
        "se_dir": DATA / "raw_output" / "niti_5nm_Ni50_single_crystal" / "se",
        "table": DATA
        / "analysis_output"
        / "niti_5nm_Ni50_single_crystal"
        / "se"
        / "se_extracted_data.txt",
    },
    "10nm": {
        "se_dir": DATA / "raw_output" / "niti_10nm_Ni50_single_crystal" / "se",
        "table": DATA
        / "analysis_output"
        / "niti_10nm_Ni50_single_crystal"
        / "se"
        / "se_extracted_data.txt",
    },
}


def find_peak_step(table_path: Path) -> int:
    arr = np.loadtxt(table_path, comments="#")
    strain = arr[:, 0]
    steps = arr[:, 8].astype(int)
    return int(steps[np.argmax(strain)])


def find_dump(se_dir: Path, step: int) -> Path:
    exact = se_dir / f"dump.se.{step}.gz"
    if exact.exists():
        return exact
    matches = sorted(se_dir.glob(f"dump.se.{step}*.gz"))
    if not matches:
        raise FileNotFoundError(f"No dump for step {step} in {se_dir}")
    return matches[0]


def counts_by_name(
    struct_array: np.ndarray, mapping: dict[int, str], mask: np.ndarray
) -> dict:
    vals, cnts = np.unique(struct_array[mask], return_counts=True)
    total = int(np.sum(mask))
    out = {}
    for v, c in zip(vals, cnts):
        name = mapping.get(int(v), f"Type-{int(v)}")
        out[name] = {
            "count": int(c),
            "fraction_within_mask": float(c / max(total, 1)),
        }
    return out


def run_case(case_name: str, cfg: dict) -> dict:
    peak_step = find_peak_step(cfg["table"])
    dump_path = find_dump(cfg["se_dir"], peak_step)

    # PTM mask defines non-B2 population used for downstream scan.
    ptm = PolyhedralTemplateMatchingModifier(rmsd_cutoff=0.15, output_orientation=False)
    p_ptm = import_file(str(dump_path))
    p_ptm.modifiers.append(ptm)
    d_ptm = p_ptm.compute()
    st_ptm = np.asarray(d_ptm.particles["Structure Type"])
    nonb2_mask = st_ptm != 3

    # CNA classification of the same frame.
    cna = CommonNeighborAnalysisModifier(
        mode=CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff
    )
    p_cna = import_file(str(dump_path))
    p_cna.modifiers.append(cna)
    d_cna = p_cna.compute()
    st_cna = np.asarray(d_cna.particles["Structure Type"])
    cna_map = {int(s.id): s.name for s in cna.structures}

    # Ackland-Jones classification of the same frame.
    aj = AcklandJonesModifier()
    p_aj = import_file(str(dump_path))
    p_aj.modifiers.append(aj)
    d_aj = p_aj.compute()
    st_aj = np.asarray(d_aj.particles["Structure Type"])
    aj_map = {int(s.id): s.name for s in aj.structures}

    return {
        "peak_step": peak_step,
        "dump_file": str(dump_path.relative_to(BASE)).replace("\\", "/"),
        "ptm_nonb2_population": int(np.sum(nonb2_mask)),
        "cna_on_ptm_nonb2": counts_by_name(st_cna, cna_map, nonb2_mask),
        "ackland_jones_on_ptm_nonb2": counts_by_name(st_aj, aj_map, nonb2_mask),
        "notes": (
            "CNA/Ackland-Jones are local-order proxies and do not uniquely separate "
            "R-phase from B19'. Results are complementary consistency checks only."
        ),
    }


def main() -> None:
    results = {name: run_case(name, cfg) for name, cfg in CASES.items()}
    out = {
        "method": (
            "Peak-frame PTM non-B2 mask followed by CNA and Ackland-Jones scans "
            "on the same atoms as additional local-order evidence."
        ),
        "cases": results,
    }
    OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(str(OUT))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
