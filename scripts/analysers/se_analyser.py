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

