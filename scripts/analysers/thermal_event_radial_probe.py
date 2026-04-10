import json
import math
from pathlib import Path

import numpy as np
from ovito.data import CutoffNeighborFinder
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier

BASE = Path(__file__).resolve().parents[2]
DATA = BASE / "data"
OUT = DATA / "analysis_output" / "uncategorized"
OUT.mkdir(parents=True, exist_ok=True)

NEIGHBOR_CUTOFF = 3.2
PTM_RMSD_CUTOFF = 0.15
CORE_FRAC = 0.80

CASES = {
    "5nm": {
        "thermal_dir": DATA
        / "raw_output"
        / "niti_5nm_Ni50_single_crystal"
        / "thermal_cycle",
        "thermal_table": DATA
        / "analysis_output"
        / "niti_5nm_Ni50_single_crystal"
        / "thermal_cycle"
        / "extracted_data.txt",
    },
    "10nm": {
        "thermal_dir": DATA
        / "raw_output"
        / "niti_10nm_Ni50_single_crystal"
        / "thermal_cycle",
        "thermal_table": DATA
        / "analysis_output"
        / "niti_10nm_Ni50_single_crystal"
        / "thermal_cycle"
        / "extracted_data.txt",
    },
}


def find_dump_step(directory: Path, step: int) -> Path:
    exact = directory / f"dump.cycle.{step}.gz"
    if exact.exists():
        return exact
    candidates = sorted(directory.glob(f"dump.cycle.{step}*.gz"))
    if not candidates:
        raise FileNotFoundError(f"No dump for step {step} in {directory}")
    return min(candidates, key=lambda p: len(p.name))


def thermal_step_near_temp(table_path: Path, target_temp: float = 272.0) -> int:
    arr = np.loadtxt(table_path, comments="#")
    temps = arr[:, 0]
    steps = arr[:, 9].astype(int)
    idx = int(np.argmin(np.abs(temps - target_temp)))
    return int(steps[idx])


def core_surface_masks(pos: np.ndarray, core_frac: float = CORE_FRAC):
    com = pos.mean(axis=0)
    d2 = np.sum((pos - com) ** 2, axis=1)
    rmax = math.sqrt(float(np.max(d2)))
    core = d2 <= (core_frac * rmax) ** 2
    return core, rmax, com


def largest_cluster(indices: np.ndarray, finder: CutoffNeighborFinder):
    idx_set = set(int(i) for i in indices.tolist())
    if not idx_set:
        return np.array([], dtype=int)

    parent = {i: i for i in idx_set}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i in idx_set:
        for neigh in finder.find(i):
            j = int(neigh.index)
            if j in idx_set and i < j:
                union(i, j)

    groups = {}
    for i in idx_set:
        r = find(i)
        groups.setdefault(r, []).append(i)

    return np.array(max(groups.values(), key=len), dtype=int)


def analyze_case(case_name: str, cfg: dict):
    step = thermal_step_near_temp(cfg["thermal_table"], 272.0)
    frame = find_dump_step(cfg["thermal_dir"], step)

    pipe = import_file(str(frame))
    pipe.modifiers.append(
        PolyhedralTemplateMatchingModifier(rmsd_cutoff=PTM_RMSD_CUTOFF)
    )
    data = pipe.compute()

    pos = np.asarray(data.particles.position)
    st = np.asarray(data.particles["Structure Type"])

    core, rmax, com = core_surface_masks(pos)
    core_nonb2 = np.flatnonzero((st != 3) & core)

    finder = CutoffNeighborFinder(NEIGHBOR_CUTOFF, data)
    cluster = largest_cluster(core_nonb2, finder)

    out = {
        "step": step,
        "frame": str(frame.relative_to(BASE)).replace("\\", "/"),
        "core_nonb2_count": int(len(core_nonb2)),
        "largest_core_nonb2_cluster": int(len(cluster)),
        "cluster_centroid_r_over_rmax": None,
    }

    if len(cluster) > 0:
        centroid = pos[cluster].mean(axis=0)
        out["cluster_centroid_r_over_rmax"] = float(
            np.linalg.norm(centroid - com) / max(rmax, 1e-12)
        )

    return out


def main():
    result = {name: analyze_case(name, cfg) for name, cfg in CASES.items()}
    out_path = OUT / "thermal_event_radial_probe.json"
    out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
