import json
import math
from pathlib import Path

import numpy as np
from ovito.data import CutoffNeighborFinder
from ovito.io import import_file
from ovito.modifiers import (
    ConstructSurfaceModifier,
    CoordinationAnalysisModifier,
    ExpressionSelectionModifier,
    PolyhedralTemplateMatchingModifier,
)

BASE = Path(__file__).resolve().parents[2]
DATA = BASE / "data"
OUT = DATA / "analysis_output" / "uncategorized"
OUT.mkdir(parents=True, exist_ok=True)

CORE_FRAC = 0.80
PTM_RMSD_CUTOFF = 0.15
NEIGHBOR_CUTOFF = 3.2
EQ_OFFSET = 1000
STEPS_LOAD = 1_000_000
END_LOAD = EQ_OFFSET + STEPS_LOAD

CASES = {
    "5nm": {
        "se_dir": DATA / "raw_output" / "niti_5nm_Ni50_single_crystal" / "se",
        "thermal_dir": DATA
        / "raw_output"
        / "niti_5nm_Ni50_single_crystal"
        / "thermal_cycle",
        "se_table": DATA
        / "analysis_output"
        / "niti_5nm_Ni50_single_crystal"
        / "se"
        / "se_extracted_data.txt",
        "thermal_table": DATA
        / "analysis_output"
        / "niti_5nm_Ni50_single_crystal"
        / "thermal_cycle"
        / "extracted_data.txt",
        "minimized": DATA
        / "raw_output"
        / "niti_5nm_Ni50_single_crystal"
        / "minimize"
        / "final_minimized.data",
    },
    "10nm": {
        "se_dir": DATA / "raw_output" / "niti_10nm_Ni50_single_crystal" / "se",
        "thermal_dir": DATA
        / "raw_output"
        / "niti_10nm_Ni50_single_crystal"
        / "thermal_cycle",
        "se_table": DATA
        / "analysis_output"
        / "niti_10nm_Ni50_single_crystal"
        / "se"
        / "se_extracted_data.txt",
        "thermal_table": DATA
        / "analysis_output"
        / "niti_10nm_Ni50_single_crystal"
        / "thermal_cycle"
        / "extracted_data.txt",
        "minimized": DATA
        / "raw_output"
        / "niti_10nm_Ni50_single_crystal"
        / "minimize"
        / "final_minimized.data",
    },
}


def load_table(path: Path) -> np.ndarray:
    return np.loadtxt(path, comments="#")


def find_dump_step(directory: Path, prefix: str, step: int) -> Path:
    exact = directory / f"{prefix}.{step}.gz"
    if exact.exists():
        return exact
    candidates = sorted(directory.glob(f"{prefix}.{step}*.gz"))
    if not candidates:
        raise FileNotFoundError(f"No dump for step {step} in {directory}")
    return min(candidates, key=lambda p: len(p.name))


def thermal_step_near_temp(thermal_table: np.ndarray, target_temp: float) -> int:
    temps = thermal_table[:, 0]
    steps = thermal_table[:, 9].astype(int)
    idx = int(np.argmin(np.abs(temps - target_temp)))
    return int(steps[idx])


def peak_and_recovered_steps(se_table: np.ndarray) -> tuple[int, int]:
    strain = se_table[:, 0]
    steps = se_table[:, 8].astype(int)
    peak_step = int(steps[np.argmax(strain)])
    post_load = se_table[se_table[:, 8] > END_LOAD]
    near_zero = post_load[np.abs(post_load[:, 0]) <= 0.002]
    if len(near_zero) == 0:
        recovered_step = int(se_table[-1, 8])
    else:
        recovered_step = int(near_zero[-1, 8])
    return peak_step, recovered_step


def build_ptm_pipeline(path: Path, with_orientation: bool = False):
    p = import_file(str(path))
    p.modifiers.append(
        PolyhedralTemplateMatchingModifier(
            rmsd_cutoff=PTM_RMSD_CUTOFF,
            output_rmsd=True,
            output_orientation=with_orientation,
        )
    )
    return p


def core_surface_masks(pos: np.ndarray, core_frac: float = CORE_FRAC):
    com = pos.mean(axis=0)
    d2 = np.sum((pos - com) ** 2, axis=1)
    rmax = math.sqrt(float(np.max(d2)))
    core = d2 <= (core_frac * rmax) ** 2
    surf = ~core
    return core, surf, rmax, com, d2


def structure_fraction_report(data, ptm_structures):
    st = np.asarray(data.particles["Structure Type"])
    n = len(st)
    out = {}
    for s in ptm_structures:
        c = int(np.count_nonzero(st == s.id))
        out[s.name] = {"count": c, "fraction": c / max(n, 1)}
    return out


def export_rdf_non_b2(case_name: str, peak_dump: Path):
    p = import_file(str(peak_dump))
    p.modifiers.append(PolyhedralTemplateMatchingModifier(rmsd_cutoff=PTM_RMSD_CUTOFF))
    p.modifiers.append(ExpressionSelectionModifier(expression="StructureType != 3"))
    p.modifiers.append(
        CoordinationAnalysisModifier(cutoff=6.0, number_of_bins=180, only_selected=True)
    )
    d = p.compute()
    tab = d.tables["coordination-rdf"]
    y = np.asarray(tab.y, dtype=float)
    n = len(y)
    start, end = tab.interval
    dr = (end - start) / n
    r = start + (np.arange(n) + 0.5) * dr

    out_csv = OUT / f"rdf_nonb2_peak_{case_name}.csv"
    with out_csv.open("w", encoding="utf-8") as f:
        f.write("r_angstrom,g_r\n")
        for rv, gv in zip(r, y):
            f.write(f"{rv:.6f},{gv:.10e}\n")

    return {
        "file": str(out_csv.relative_to(BASE)).replace("\\", "/"),
        "cutoff_angstrom": 6.0,
        "bins": int(n),
        "peak_g_r": float(np.max(y)),
        "peak_r_angstrom": float(r[int(np.argmax(y))]),
    }


def nonb2_angle_histogram(data, bins: int = 180, cutoff: float = NEIGHBOR_CUTOFF):
    st = np.asarray(data.particles["Structure Type"])
    pos = np.asarray(data.particles.position)
    selected = np.flatnonzero(st != 3)

    finder = CutoffNeighborFinder(cutoff, data)
    angles = []

    for i in selected:
        vecs = []
        for neigh in finder.find(i):
            vecs.append(np.array([neigh.delta[0], neigh.delta[1], neigh.delta[2]]))
        if len(vecs) < 2:
            continue
        for a in range(len(vecs) - 1):
            va = vecs[a]
            na = np.linalg.norm(va)
            if na == 0:
                continue
            for b in range(a + 1, len(vecs)):
                vb = vecs[b]
                nb = np.linalg.norm(vb)
                if nb == 0:
                    continue
                cosang = float(np.dot(va, vb) / (na * nb))
                cosang = max(-1.0, min(1.0, cosang))
                angles.append(math.degrees(math.acos(cosang)))

    if not angles:
        hist = np.zeros(bins, dtype=float)
        centers = np.linspace(0.5, 179.5, bins)
        return centers, hist, 0

    hist, edges = np.histogram(angles, bins=bins, range=(0.0, 180.0), density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, hist, len(angles)


def export_nonb2_bond_angles(case_name: str, peak_dump: Path):
    p = build_ptm_pipeline(peak_dump)
    d = p.compute()
    centers, hist, n_angles = nonb2_angle_histogram(d)

    out_csv = OUT / f"bond_angle_nonb2_peak_{case_name}.csv"
    with out_csv.open("w", encoding="utf-8") as f:
        f.write("theta_deg,pdf\n")
        for t, h in zip(centers, hist):
            f.write(f"{t:.6f},{h:.10e}\n")

    peak_indices = np.argsort(hist)[-3:][::-1]
    peak_angles = [float(centers[i]) for i in peak_indices]

    return {
        "file": str(out_csv.relative_to(BASE)).replace("\\", "/"),
        "angle_samples": int(n_angles),
        "top_peak_angles_deg": peak_angles,
        "neighbor_cutoff_angstrom": NEIGHBOR_CUTOFF,
    }


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

    buckets = {}
    for i in idx_set:
        r = find(i)
        buckets.setdefault(r, []).append(i)

    best = max(buckets.values(), key=len)
    return np.array(best, dtype=int)


def nucleation_tracking(case_name: str, se_dir: Path, se_table: np.ndarray):
    rows = se_table[(se_table[:, 8] >= EQ_OFFSET) & (se_table[:, 8] <= END_LOAD)]
    rows = rows[np.argsort(rows[:, 8])]

    track = []
    first_event = None

    for row in rows:
        strain = float(row[0])
        step = int(row[8])
        dump = find_dump_step(se_dir, "dump.se", step)

        p = build_ptm_pipeline(dump)
        d = p.compute()

        pos = np.asarray(d.particles.position)
        st = np.asarray(d.particles["Structure Type"])
        core, _, rmax, com, d2 = core_surface_masks(pos)

        core_nonb2 = np.flatnonzero((st != 3) & core)
        finder = CutoffNeighborFinder(NEIGHBOR_CUTOFF, d)
        cluster = largest_cluster(core_nonb2, finder)

        if len(cluster) > 0:
            centroid = pos[cluster].mean(axis=0)
            r_ratio = float(np.linalg.norm(centroid - com) / max(rmax, 1e-12))
        else:
            r_ratio = None

        track.append(
            {
                "step": step,
                "strain": strain,
                "core_nonb2_count": int(len(core_nonb2)),
                "largest_core_nonb2_cluster": int(len(cluster)),
                "cluster_centroid_r_over_rmax": r_ratio,
            }
        )

        if first_event is None and len(cluster) >= 20:
            first_event = track[-1]

    out_csv = OUT / f"nucleation_tracking_{case_name}.csv"
    with out_csv.open("w", encoding="utf-8") as f:
        f.write(
            "step,strain,core_nonb2_count,largest_core_nonb2_cluster,cluster_centroid_r_over_rmax\n"
        )
        for r in track:
            rr = (
                ""
                if r["cluster_centroid_r_over_rmax"] is None
                else f"{r['cluster_centroid_r_over_rmax']:.6f}"
            )
            f.write(
                f"{r['step']},{r['strain']:.8f},{r['core_nonb2_count']},{r['largest_core_nonb2_cluster']},{rr}\n"
            )

    return {
        "file": str(out_csv.relative_to(BASE)).replace("\\", "/"),
        "first_cluster_event_threshold_atoms": 20,
        "first_event": first_event,
    }


def rmsd_core_surface_equilibrium(se_dir: Path):
    dump = find_dump_step(se_dir, "dump.se", 0)
    p = build_ptm_pipeline(dump)
    d = p.compute()
    pos = np.asarray(d.particles.position)
    core, surf, _, _, _ = core_surface_masks(pos)

    if "RMSD" not in d.particles:
        return {"error": "PTM RMSD property not available"}

    rmsd = np.asarray(d.particles["RMSD"], dtype=float)

    def stats(vals):
        if len(vals) == 0:
            return {}
        return {
            "n": int(len(vals)),
            "mean": float(np.mean(vals)),
            "median": float(np.median(vals)),
            "p90": float(np.percentile(vals, 90)),
            "p99": float(np.percentile(vals, 99)),
        }

    return {
        "core": stats(rmsd[core]),
        "surface": stats(rmsd[surf]),
        "rmsd_cutoff": PTM_RMSD_CUTOFF,
    }


def find_stress_zz(data):
    keys = list(data.particles.keys())
    if "c_st[3]" in data.particles:
        arr = np.asarray(data.particles["c_st[3]"], dtype=float)
        return arr
    if "c_st" in data.particles:
        st = np.asarray(data.particles["c_st"], dtype=float)
        if st.ndim == 2 and st.shape[1] >= 3:
            return st[:, 2]
    for k in keys:
        sk = str(k).lower()
        if "c_st" in sk and ("[3]" in sk or sk.endswith("3")):
            arr = np.asarray(data.particles[k], dtype=float)
            if arr.ndim == 1:
                return arr
    return None


def stress_volume_sensitivity(case_name: str, se_dir: Path, se_table: np.ndarray):
    top = se_table[np.argsort(se_table[:, 1])[-5:]]
    top = top[np.argsort(top[:, 8])]

    rows = []
    for r in top:
        step = int(r[8])
        dump = find_dump_step(se_dir, "dump.se", step)

        p = import_file(str(dump))
        p.modifiers.append(
            ConstructSurfaceModifier(
                radius=2.5, smoothing_level=0, identify_regions=True
            )
        )
        d = p.compute()

        pos = np.asarray(d.particles.position)
        _, _, rmax, _, _ = core_surface_masks(pos)
        v_rmax = (4.0 / 3.0) * math.pi * (rmax**3)
        v_mesh = float(d.attributes.get("ConstructSurfaceMesh.filled_volume", np.nan))

        szz = find_stress_zz(d)
        if szz is None:
            sigma_rmax = float("nan")
            sigma_mesh = float("nan")
        else:
            sigma_rmax = -float(np.sum(szz)) / max(v_rmax, 1e-12) * 1.0e-4
            sigma_mesh = -float(np.sum(szz)) / max(v_mesh, 1e-12) * 1.0e-4

        rel = float("nan")
        if (
            np.isfinite(sigma_rmax)
            and abs(sigma_rmax) > 1e-12
            and np.isfinite(sigma_mesh)
        ):
            rel = 100.0 * (sigma_mesh - sigma_rmax) / sigma_rmax

        rows.append(
            {
                "step": step,
                "sigma_rmax_gpa": sigma_rmax,
                "sigma_mesh_gpa": sigma_mesh,
                "mesh_vs_rmax_percent": rel,
            }
        )

    vals = [
        x["mesh_vs_rmax_percent"]
        for x in rows
        if np.isfinite(x["mesh_vs_rmax_percent"])
    ]
    summary = {
        "n_frames": len(rows),
        "mean_percent": float(np.mean(vals)) if vals else float("nan"),
        "max_abs_percent": float(np.max(np.abs(vals))) if vals else float("nan"),
    }

    out_csv = OUT / f"stress_norm_sensitivity_{case_name}.csv"
    with out_csv.open("w", encoding="utf-8") as f:
        f.write("step,sigma_rmax_gpa,sigma_mesh_gpa,mesh_vs_rmax_percent\n")
        for rr in rows:
            f.write(
                f"{rr['step']},{rr['sigma_rmax_gpa']:.8f},{rr['sigma_mesh_gpa']:.8f},{rr['mesh_vs_rmax_percent']:.6f}\n"
            )

    return {
        "file": str(out_csv.relative_to(BASE)).replace("\\", "/"),
        "summary": summary,
        "surface_mesh_radius_angstrom": 2.5,
    }


def loop_area_metrics(se_table: np.ndarray):
    load = se_table[(se_table[:, 8] >= EQ_OFFSET) & (se_table[:, 8] <= END_LOAD)]
    unload = se_table[se_table[:, 8] > END_LOAD]
    load = load[np.argsort(load[:, 0])]
    unload = unload[np.argsort(unload[:, 0])]

    eps_load = load[:, 0]
    sig_load = load[:, 1]

    eps_un = unload[:, 0][::-1]
    sig_un = unload[:, 1][::-1]

    w_load = float(np.trapezoid(sig_load, eps_load))
    w_unload = float(np.trapezoid(sig_un, eps_un))

    eps_common = np.linspace(0.0, min(eps_load.max(), eps_un.max()), 500)
    sig_load_i = np.interp(eps_common, eps_load, sig_load)
    sig_un_i = np.interp(eps_common, eps_un, sig_un)
    loop = float(np.trapezoid(sig_load_i - sig_un_i, eps_common))

    peak_stress = float(np.max(sig_load))
    peak_strain = float(np.max(eps_load))
    norm = float(loop / max(peak_stress * peak_strain, 1e-12))

    return {
        "w_loading_gpa": w_load,
        "w_unloading_gpa": w_unload,
        "loop_area_gpa": loop,
        "normalized_loop_area": norm,
        "peak_stress_gpa": peak_stress,
        "peak_strain": peak_strain,
    }


def minimization_ptm_check(min_path: Path):
    p = build_ptm_pipeline(min_path, with_orientation=True)
    d = p.compute()
    st = np.asarray(d.particles["Structure Type"])

    bcc_frac = float(np.count_nonzero(st == 3) / max(len(st), 1))

    result = {"bcc_fraction": bcc_frac}
    if "Orientation" in d.particles:
        ori = np.asarray(d.particles["Orientation"])
        result["orientation_quaternion_component_std"] = [
            float(np.std(ori[:, i])) for i in range(min(4, ori.shape[1]))
        ]
    return result


def main():
    report = {
        "ptm_r_phase_support": None,
        "ptm_keyframe_phase_fractions": {},
        "rdf_nonb2_peak": {},
        "bond_angle_nonb2_peak": {},
        "nucleation_tracking": {},
        "ptm_rmsd_core_surface_equilibrium": {},
        "stress_norm_sensitivity": {},
        "loop_area_metrics": {},
        "minimization_ptm_check": {},
        "notes": [],
    }

    ptm_probe = PolyhedralTemplateMatchingModifier(rmsd_cutoff=PTM_RMSD_CUTOFF)
    names = [s.name.lower() for s in ptm_probe.structures]
    report["ptm_r_phase_support"] = {
        "supported": any("r-phase" in n or "r phase" in n or n == "r" for n in names),
        "available_templates": [s.name for s in ptm_probe.structures],
    }

    if not report["ptm_r_phase_support"]["supported"]:
        report["notes"].append(
            "Current OVITO PTM build has no explicit R-phase template; R-phase cannot be directly quantified via PTM Structure Type."
        )

    for case_name, cfg in CASES.items():
        se_tab = load_table(cfg["se_table"])
        th_tab = load_table(cfg["thermal_table"])

        peak_step, recovered_step = peak_and_recovered_steps(se_tab)
        t272_step = thermal_step_near_temp(th_tab, 272.0)

        frames = {
            "equilibrium_500K": find_dump_step(cfg["thermal_dir"], "dump.cycle", 0),
            "cooling_near_272K": find_dump_step(
                cfg["thermal_dir"], "dump.cycle", t272_step
            ),
            "peak_strain": find_dump_step(cfg["se_dir"], "dump.se", peak_step),
            "recovered": find_dump_step(cfg["se_dir"], "dump.se", recovered_step),
        }

        case_ptm = {}
        for label, frame in frames.items():
            p = build_ptm_pipeline(frame)
            d = p.compute()
            case_ptm[label] = {
                "step": int(frame.stem.split(".")[-1].split("(")[0]),
                "file": str(frame.relative_to(BASE)).replace("\\", "/"),
                "fractions": structure_fraction_report(d, p.modifiers[0].structures),
            }
        report["ptm_keyframe_phase_fractions"][case_name] = case_ptm

        report["rdf_nonb2_peak"][case_name] = export_rdf_non_b2(
            case_name, frames["peak_strain"]
        )
        report["bond_angle_nonb2_peak"][case_name] = export_nonb2_bond_angles(
            case_name, frames["peak_strain"]
        )
        report["nucleation_tracking"][case_name] = nucleation_tracking(
            case_name, cfg["se_dir"], se_tab
        )
        report["ptm_rmsd_core_surface_equilibrium"][case_name] = (
            rmsd_core_surface_equilibrium(cfg["se_dir"])
        )
        report["stress_norm_sensitivity"][case_name] = stress_volume_sensitivity(
            case_name, cfg["se_dir"], se_tab
        )
        report["loop_area_metrics"][case_name] = loop_area_metrics(se_tab)
        report["minimization_ptm_check"][case_name] = minimization_ptm_check(
            cfg["minimized"]
        )

    out_json = OUT / "revision_prewrite_analysis.json"
    out_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(str(out_json))


if __name__ == "__main__":
    main()
