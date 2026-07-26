                      
"""
OVITO Python Rendering Script for NiTi Nanoparticle MD Paper
=============================================================
Generates atomistic renderings using OVITO's Python API + Tachyon renderer.
Requires: ovito (Python package), pip install ovito

Run from repo root:  python scripts/paper/generate_ovito_renderings.py

Outputs to: paper/figures/

Figure Plan:
  Fig 1(a,b) - Equilibrated NPs (PTM colored, whole + cross-section), 5 nm
  Fig 1(c,d) - Equilibrated NPs (PTM colored, whole + cross-section), 10 nm
  Fig 7      - SE snapshots at key strains: ε=0%, 5%, 10%, 5%(unload), 0%(recovered)

"""

import os
import sys
import glob
import re
import math
import traceback
import warnings

import numpy as np

try:
                                                                                
                                                                              
                                          
    warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")

                                                      
    import ovito._extensions.pyscript              

    from ovito.io import import_file
    from ovito.modifiers import (
        PolyhedralTemplateMatchingModifier,
        SliceModifier,
        ExpressionSelectionModifier,
        ColorCodingModifier,
        PythonScriptModifier,
    )
    from ovito.vis import (
        Viewport,
        TachyonRenderer,
        ParticlesVis,
    )
    from ovito.data import DataCollection
except Exception as exc:
    print("ERROR: Failed to import OVITO's Python module.")
    print("\nRoot cause (full traceback):")
    traceback.print_exc()
    print(
        "\nMost common fix on Windows: don't use the PyPI 'ovito' package inside an Anaconda base env.\n"
        "Instead use ONE of these supported setups:\n\n"
        "A) New conda env + OVITO conda package (recommended):\n"
        "  conda create -n ovito310 python=3.10 -y\n"
        "  conda activate ovito310\n"
        "  conda install -c ovito ovito -y\n"
        "  python scripts/paper/generate_ovito_renderings.py\n\n"
        "B) python.org + venv (no Anaconda):\n"
        "  py -3.10 -m venv .venv-ovito\n"
        "  .\\.venv-ovito\\Scripts\\activate\n"
        "  pip install ovito\n"
        "  python scripts/paper/generate_ovito_renderings.py\n\n"
        "If you want me to wire a VS Code interpreter + task for option (A), tell me which env name you prefer."
    )
    sys.exit(1)

                    
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUTDIR = os.path.join(REPO_ROOT, "paper", "figures")
os.makedirs(OUTDIR, exist_ok=True)

RAW_BASE = os.path.join(REPO_ROOT, "data", "raw_output")

                                                                           
EQUILIBRIUM_DUMPS = {
    "5nm": os.path.join(RAW_BASE, "niti_5nm_Ni50_single_crystal", "se", "dump.se.0.gz"),
    "10nm": os.path.join(
        RAW_BASE, "niti_10nm_Ni50_single_crystal", "se", "dump.se.0.gz"
    ),
}

                     
SE_DIRS = {
    "5nm": os.path.join(RAW_BASE, "niti_5nm_Ni50_single_crystal", "se"),
    "10nm": os.path.join(RAW_BASE, "niti_10nm_Ni50_single_crystal", "se"),
}

                                   
RMSD = 0.15
STRAIN_MAX = 0.10
STEPS_LOAD = 1000000
EQ_OFFSET = 1000

                                 
                            
                                                                      
                                             
PTM_COLORS = {
    0: (0.75, 0.75, 0.75, 1.0),                          
    1: (0.12, 0.47, 0.71, 1.0),                    
    2: (1.00, 0.50, 0.05, 1.0),                      
    3: (0.17, 0.63, 0.17, 1.0),                     
    4: (0.84, 0.15, 0.16, 1.0),                   
    5: (0.58, 0.40, 0.74, 1.0),                      
}

PARTICLE_RADIUS = 0.8                   


def step_from_strain(strain):
    """Convert target strain to nearest dump timestep."""
    step = int(EQ_OFFSET + (strain / STRAIN_MAX) * STEPS_LOAD)
                                           
    step = round(step / 5000) * 5000
    return step


def step_unload_from_strain(strain):
    """Convert target strain during unloading to nearest dump timestep."""
    END_LOAD = EQ_OFFSET + STEPS_LOAD
    frac_remaining = strain / STRAIN_MAX
    step = int(END_LOAD + (1.0 - frac_remaining) * STEPS_LOAD)
    step = round(step / 5000) * 5000
    return step


def find_dump_file(directory, step):
    """Find the dump file for a given step, handling duplicates like (1)."""
    exact = os.path.join(directory, f"dump.se.{step}.gz")
    if os.path.exists(exact):
        return exact
                                  
    files = glob.glob(os.path.join(directory, "dump.se.*.gz"))
    best = None
    best_diff = float("inf")
    for f in files:
        m = re.search(r"dump\.se\.(\d+)(?:\(\d+\))?\.gz", os.path.basename(f))
        if m:
            s = int(m.group(1))
            if abs(s - step) < best_diff:
                best_diff = abs(s - step)
                best = f
    return best


def apply_ptm_colors(frame, data):
    """Modifier function: set particle colors based on PTM Structure Type.

    This is used as a PythonScriptModifier so colors are applied in the
    pipeline (and therefore visible during rendering), not on a detached copy.
    """
    structs = np.asarray(data.particles["Structure Type"])
    n = len(structs)
    colors = np.zeros((n, 3), dtype=float)
    for stype, rgba in PTM_COLORS.items():
        mask = structs == stype
        colors[mask] = rgba[:3]

    data.particles_.create_property("Color", data=colors)


def setup_pipeline(dump_path, slice_np=False):
    """Create an OVITO pipeline with PTM and optional cross-section slice."""
    pipeline = import_file(dump_path)

                                                                
    if pipeline.source.data.cell is not None:
        pipeline.source.data.cell.vis.enabled = False

                  
    ptm = PolyhedralTemplateMatchingModifier(rmsd_cutoff=RMSD, output_orientation=False)
    pipeline.modifiers.append(ptm)

                                                                          
    pipeline.modifiers.append(PythonScriptModifier(function=apply_ptm_colors))

    if slice_np:
                                                                       
                                                                           
        data_tmp = pipeline.compute()
        pos = np.asarray(data_tmp.particles.position)
        com_y = float(np.mean(pos[:, 1]))

                                                                   
        r_from_com = np.sqrt(np.sum((pos - np.mean(pos, axis=0)) ** 2, axis=1))
        diameter = 2.0 * np.max(r_from_com)
        slab_w = max(10.0, 0.20 * diameter)                 

        slicer = SliceModifier(
            distance=com_y,
            normal=(0, 1, 0),               
            slab_width=slab_w,
        )
        pipeline.modifiers.append(slicer)

    return pipeline


def render_viewport(
    pipeline, output_path, size_px=(1600, 1600), camera_dir=(0, 0, -1), zoom_factor=1.0
):
    """Render the pipeline from a given camera direction."""
                                                          
    vis = pipeline.source.data.particles.vis
    if vis is not None:
        vis.radius = PARTICLE_RADIUS
        vis.shape = ParticlesVis.Shape.Sphere

                                                                              
                                                             
    pipeline.add_to_scene()
    try:
                                                       
        pipeline.compute()

                         
        vp = Viewport(type=Viewport.Type.Ortho, camera_dir=camera_dir)
        vp.zoom_all(size=size_px)

                    
        if zoom_factor != 1.0:
            vp.fov *= zoom_factor

                             
        renderer = TachyonRenderer(
            shadows=False,
            ambient_occlusion=True,
            ambient_occlusion_brightness=0.8,
            antialiasing_samples=12,
            direct_light_intensity=0.9,
        )

        vp.render_image(
            size=size_px,
            filename=output_path,
            renderer=renderer,
            background=(1, 1, 1),                    
        )
        print(f"  Rendered: {os.path.basename(output_path)}")
    finally:
        pipeline.remove_from_scene()


                            
def fig1_equilibrium_renderings():
    """Render equilibrated NPs: whole view and cross-section for both sizes."""
    print("\n--- Figure 1: Equilibrated nanoparticle renderings ---")

    for size in ["5nm", "10nm"]:
        dump_path = EQUILIBRIUM_DUMPS[size]
        if not os.path.exists(dump_path):
            print(f"  WARNING: Missing {dump_path}, skipping {size}")
            continue

                                  
        pipeline = setup_pipeline(dump_path, slice_np=False)
        render_viewport(
            pipeline,
            os.path.join(OUTDIR, f"fig1_{size}_whole.png"),
            camera_dir=(2, 1, -3),
            zoom_factor=1.1,
        )

                             
        pipeline = setup_pipeline(dump_path, slice_np=True)
        render_viewport(
            pipeline,
            os.path.join(OUTDIR, f"fig1_{size}_cross_section.png"),
            camera_dir=(0, -1, 0),
            zoom_factor=1.1,
        )


                                       
def fig7_se_snapshots():
    """Render cross-section snapshots at key strains during SE loading/unloading."""
    print("\n--- Figure 7: SE cross-section snapshots ---")

                                                                                 
    targets = [
        ("eq", 0.00, "load"),
        ("5pct_load", 0.05, "load"),
        ("10pct_peak", 0.10, "load"),
        ("5pct_unload", 0.05, "unload"),
        ("0pct_recovered", 0.00, "unload"),
    ]

    for size in ["5nm", "10nm"]:
        se_dir = SE_DIRS[size]

        for label, strain, phase in targets:
            if phase == "load":
                if strain == 0.0:
                    step = 0
                else:
                    step = step_from_strain(strain)
            else:
                if strain == 0.0:
                    step = 2000000              
                else:
                    step = step_unload_from_strain(strain)

            dump_path = find_dump_file(se_dir, step)
            if dump_path is None:
                print(f"  WARNING: No dump for {size} {label} (step ~{step})")
                continue

            pipeline = setup_pipeline(dump_path, slice_np=True)
            render_viewport(
                pipeline,
                os.path.join(OUTDIR, f"fig7_{size}_{label}.png"),
                camera_dir=(0, -1, 0),
                zoom_factor=1.1,
            )


                        
def generate_color_legend():
    """Create a standalone color legend image for PTM structure types."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    fig, ax = plt.subplots(figsize=(3, 1.5))
    ax.set_axis_off()

    legend_items = [
        ("B2 (austenite)", PTM_COLORS[3][:3]),
        ("Other", PTM_COLORS[0][:3]),
        ("FCC", PTM_COLORS[1][:3]),
        ("HCP", PTM_COLORS[2][:3]),
        ("ICO", PTM_COLORS[4][:3]),
    ]

    patches = [mpatches.Patch(color=c, label=l) for l, c in legend_items]
    ax.legend(
        handles=patches,
        loc="center",
        ncol=2,
        frameon=False,
        fontsize=9,
        handlelength=1.5,
        handleheight=1.2,
    )

    fig.savefig(
        os.path.join(OUTDIR, "ptm_color_legend.pdf"),
        bbox_inches="tight",
        pad_inches=0.05,
        dpi=300,
    )
    fig.savefig(
        os.path.join(OUTDIR, "ptm_color_legend.png"),
        bbox_inches="tight",
        pad_inches=0.05,
        dpi=300,
    )
    plt.close(fig)
    print("  [OK] PTM color legend")


      
if __name__ == "__main__":
    print("=" * 60)
    print("OVITO Python Rendering Pipeline")
    print(f"Output directory: {OUTDIR}")
    print("=" * 60)

    fig1_equilibrium_renderings()
    fig7_se_snapshots()
    generate_color_legend()

    print("=" * 60)
    print("All OVITO renderings complete.")
    print("=" * 60)
