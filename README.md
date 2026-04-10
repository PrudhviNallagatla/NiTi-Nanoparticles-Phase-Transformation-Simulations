# NiTi Nanoparticle SME/SE Simulation Workspace

This repository houses a high-throughput computational workflow for simulating **Nickel-Titanium (NiTi)** nanoparticles using Molecular Dynamics (MD). It is designed to investigate two critical material behaviors—the **Shape Memory Effect (SME)** and **Superelasticity (SE)**—under realistic synthesis conditions.

## 📂 Repository Structure

The workspace keeps the simulation inputs, analysis scripts, manuscript assets, and generated outputs separate. Large dump, restart, and analysis-result files are intentionally collapsed here.

```text
.
├── all_dependencies.sh
├── log.lammps
├── README.md
├── docs/
│   ├── acta_revision_master_todo.md
│   ├── Azure_workflow.md
│   ├── NiTi_Nanoparticle_MD_Report.md
│   ├── response_letter_draft.md
│   └── reviewer_response_matrix.md
├── paper/
│   ├── main.tex
│   ├── paper_for_reference.tex
│   ├── references.bib
│   ├── elsarticle.cls
│   ├── elsarticle-harv.bst
│   ├── elsarticle-num.bst
│   ├── elsarticle-num-names.bst
│   └── figures/
├── potentials/
│   ├── Ni-Ti.eam.fs
│   ├── NiTi.meam
│   ├── library.meam
│   ├── NiTi_graph_v1.pb
│   └── NiTi_graph_v2.pb
├── scripts/
│   ├── run_size_predictor.sh
│   ├── se.sh
│   ├── sme.sh
│   ├── sme_and_se.sh
│   ├── structure_gen_ase.sh
│   ├── thermal_cycle_script.sh
│   ├── zip_it_all.sh
│   ├── analysers/
│   │   ├── se_analyser.py
│   │   ├── thermal_analyser.py
│   │   ├── thermal_event_radial_probe.py
│   │   ├── revision_prewrite_analysis.py
│   │   ├── ovito_revision_checks.py
│   │   ├── rphase_proxy_analysis.py
│   │   ├── rphase_cna_aj_scan.py
│   │   └── .legacy_scripts/
│   │       ├── hysteresis_analyser.py
│   │       └── analyze_thermal_cycle.py
│   └── paper/
│       ├── generate_data_figures.py
│       ├── generate_high_value_figures.py
│       ├── generate_ovito_renderings.py
│       ├── generate_radial_profile.py
│       └── generate_revision_figures.py
├── src/
│   ├── ase_np_gen.py
│   ├── niti_np_gen.py
│   └── lmps/
│       ├── minimize.lmp
│       ├── se.lmp
│       ├── sme.lmp
│       └── thermal_cycle.lmp
├── data/
│   ├── structures/
│   │   ├── niti_5nm_Ni50_single_crystal.xyz
│   │   └── niti_10nm_Ni50_single_crystal.xyz
│   ├── raw_output/
│   │   ├── niti_2nm_Ni50_amorphous/
│   │   ├── niti_5nm_Ni50_single_crystal/
│   │   └── niti_10nm_Ni50_single_crystal/
│   └── analysis_output/
│       ├── niti_2nm_Ni50_amorphous/
│       ├── niti_5nm_Ni50_single_crystal/
│       ├── niti_10nm_Ni50_single_crystal/
│       └── uncategorized/
├── logs/
│   ├── 2nm_1grain/
│   └── 5nm_1grain/
├── zips/
└── .devcontainer/
```

* **`src/lmps/`**: LAMMPS input scripts for minimization, thermal cycling, SME, and SE.
* **`scripts/analysers/`**: Python post-processing scripts for phase analysis, thermal analysis, and paper figures.
* **`scripts/paper/`**: Figure-generation scripts used to build the manuscript graphics.
* **`paper/`**: LaTeX manuscript, bibliography, and class/style files.
* **`data/`**: Input structures plus raw and processed simulation outputs.
* **`potentials/`**: Interatomic potential files and related graph models.
* **`docs/`**: Revision notes, response materials, and manuscript support documents.

Generated dump, restart, and analysis files are omitted from the tree above for readability; the directory structure is preserved where those artifacts live.

---

## ⚙️ Simulation Protocol (How to Proceed)

This workflow is sequential. You cannot run the mechanical tests (SME/SE) until you have successfully established the thermal baseline.

### Step 1: Establish the Thermal Baseline (`thermal_cycle`)

**What it does:**
This simulation takes your generated nanoparticle and runs it through a full heating and cooling loop. It heats the particle to a high temperature (Austenite phase) and cools it to a low temperature (Martensite phase).

**Why it is needed:**
Real materials have a "thermal history." You cannot simply start a simulation at 0 Kelvin. This step generates "Restart Files"—snapshots of the atoms at critical points (fully cooled, fully heated)—which serve as the starting points for the next steps.

**Action:**
Run the thermal cycle simulation first. Wait for it to complete the full loop (Heating  Cooling  Heating).

### Step 2: Validate the Physics (Critical Analysis)

**What it is:**
Before running mechanical tests, you must verify that your particle actually transformed.

* ** (Martensite Finish):** The temperature where the particle becomes fully Martensite.
* ** (Austenite Finish):** The temperature where the particle becomes fully Austenite.

**How to do it:**
Use the analysis tools to look at the "Phase Fraction" from Step 1.

1. Check if the particle is **100% Martensite** at your target low temperature (e.g., 200K).
2. Check if the particle is **100% Austenite** at your target high temperature (e.g., 500K).

**Decision:**

* If the particle is *not* fully Martensite at 200K, you must lower the temperature in the SME script.
* If the particle is *not* fully Austenite at 500K, you must raise the temperature in the SE script.

### Step 3: Mechanical Testing

Once the thermal baseline is verified, you can run the following two simulations. They are independent and can be run simultaneously.

#### A. Shape Memory Effect (`sme`)

**The Goal:** Demonstrate that the material can "remember" its original shape.
**The Process:**

1. **Load:** The simulation picks up the "Fully Cooled" snapshot from Step 1.
2. **Deform:** It mechanically crushes (compresses) the particle at low temperature. The material stays deformed after unloading.
3. **Recover:** It heats the particle back up. If successful, the particle will spontaneously snap back to its original shape.

#### B. Superelasticity (`se`)

**The Goal:** Demonstrate rubber-like elasticity at high temperatures.
**The Process:**

1. **Load:** The simulation picks up the "Fully Heated" snapshot from Step 1.
2. **Cycle:** It compresses and releases the particle at high temperature.
3. **Result:** unlike normal metals, it should recover massive strains immediately upon unloading, producing a characteristic "flag-shaped" hysteresis loop in the data.

---

## 💻 Hardware & Execution Environment

These simulations are computationally expensive due to the complexity of the MEAM potential. A standard laptop is generally insufficient for production runs.

### Recommended Hardware

* **CPU Architecture:**
* **Minimum:** 8 Physical Cores (e.g., Modern Ryzen 7 / Intel i7).
* **Recommended:** 16-32 Physical Cores (Workstation or HPC Node).
* *Note:* The MEAM potential scales reasonably well with MPI, but performance gains diminish beyond ~64 cores for small nanoparticles (<50k atoms).


* **Memory (RAM):**
* **Minimum:** 16 GB.
* **Recommended:** 32 GB or higher. Large atomic dump files and trajectory analysis require significant memory overhead.


* **Storage:**
* **Requirement:** High-speed SSD.
* **Capacity:** Ensure at least **50-100 GB** of free space. A full SME cycle generates thousands of frame snapshots; insufficient disk space will cause the simulation to crash mid-run.



### Software Prerequisites

To run these scripts, your environment must have the following installed and configured:

1. **LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator)**
* **Version:** Stable release (2020 or later).
* **Required Packages:** You *must* compile LAMMPS with the following packages enabled:
* `MEAM` (Modified Embedded Atom Method)
* `MANYBODY` (For multi-body interactions)
* `MPI` (Message Passing Interface for parallel execution)
* `MOLECULE` (Standard atomistic definitions)




2. **Python Environment**
* **Python 3.8+**
* **Libraries:** `numpy` (math), `pandas` (data handling), `scipy` (curve fitting).



---

## 📚 References

This project relies on validated interatomic potentials and methodologies from the following literature:

* **Primary Potential (Ni-Ti MEAM):**
* Ko, W. S., Grabowski, B., & Neugebauer, J. (2015). Development of a modified embedded-atom method potential for Ni-Ti systems. *Physical Review B*, 92(13), 134107.


* **Structure Identification:**
* Ackland, G. J., & Jones, A. P. (2006). Applications of local crystal structure measures in atomistic simulations. *Physical Review B*, 73(5), 054104.


* **Phase Transformation Methodology:**
* Waitz, T., et al. (2004). Size-dependent martensitic transformation in NiTi nanoparticles. *Acta Materialia*, 52(1), 137-147.

