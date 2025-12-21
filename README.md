# NiTi Nanoparticle SME/SE Simulation Workspace

This repository houses a high-throughput computational workflow for simulating **Nickel-Titanium (NiTi)** nanoparticles. It is designed to investigate the **Shape Memory Effect (SME)** and **Superelasticity (SE)** using Molecular Dynamics (MD).

The workspace specifically targets the simulation of realistic **Electrical Discharge Machining (EDM)** synthesis conditions, enabling the study of rapid quenching effects, amorphous structures, and size-dependent phase transformations.

---

## 📂 Repository Architecture

The workspace follows a strict separation of concerns, decoupling the simulation logic from the physical models and the generated data. This modular structure ensures scalability and data integrity.

```text
.
├── scripts/                      # ORCHESTRATION LAYER (The Controller)
│   ├── run_all.sh                # Master pipeline trigger
│   ├── structure_gen.sh          # Batch generator wrapper
│   ├── se_pipeline.sh            # Superelasticity workflow logic
│   └── sme_scripts/              # Shape Memory Effect workflow logic
│
├── src/                          # SOURCE CODE LAYER (The Engine)
│   ├── niti_np_gen.py            # Python-based nanoparticle generator
│   └── lmps/                     # LAMMPS input script templates
│
├── data/                         # DATA STORAGE LAYER (The Vault)
│   ├── structures/               # Generated input files (.data)
│   ├── raw_output/               # Simulation dumps and thermodynamic logs
│   └── analysis_output/          # Final plots and calculated properties
│
├── potentials/                   # PHYSICS LAYER (The Model)
│   ├── NiTi.meam                 # MEAM parameter file
│   └── library.meam              # Element library
│
└── logs/                         # LOGGING

```

### Component Breakdown

* **Orchestration Layer (`scripts/`)**: Shell scripts that drive the entire lifecycle. These scripts manage file paths, handle parallel execution, and trigger the physics engine.
* **Source Layer (`src/`)**: Contains the procedural generation logic for creating nanoparticles and the raw LAMMPS templates.
* **Physics Layer (`potentials/`)**: Isolates the interatomic potentials to ensure modularity.
* **Data Layer (`data/`)**: A centralized hub for all simulation I/O, ensuring that raw data and analyzed results remain distinct from the codebase.

---

## ⚙️ Workflow Overview

The pipeline operates in three sequential stages:

1. **Structure Generation**: The Python generator (`src/`) constructs the nanoparticle geometry, applying specific EDM characteristics such as amorphous fractions or surface defects.
2. **Simulation Execution**: The orchestration scripts (`scripts/`) detect new structures and automatically queue them for processing. This includes:
* **SME Pipeline**: Thermal cycling to determine transformation temperatures ().
* **SE Pipeline**: Mechanical loading cycles to observe the superelastic stress-strain hysteresis.


3. **Analysis**: Post-processing notebooks utilize the raw thermodynamic logs to extract material properties and generate visualizations.

---

## 💻 System Requirements

This workflow is computationally intensive. The following specifications are recommended to ensure reasonable simulation times and data handling capabilities.

### Hardware Recommendations

* **Processor (CPU)**: Minimum **8 Physical Cores** recommended.
* *Note*: LAMMPS scales well with MPI; running on fewer cores will significantly increase runtime.


* **Memory (RAM)**: Minimum **16 GB**.
* *Reasoning*: Large-scale atom dumps and post-processing of thermodynamic logs require substantial memory overhead.


* **Storage**: SSD with high I/O throughput.
* *Space*: Ensure at least 50GB of free space for raw dump files before data reduction.



### Software Prerequisites

* **LAMMPS**: Compiled with `MPI` support and the `MEAM` package installed.
* **Python 3.8+**: Required for the generation and analysis modules.
* **Python Libraries**: `numpy`, `pandas`, `scipy`, `matplotlib`, `jupyter`.

<!-- ---

## ⏱️ Performance Benchmarks

Simulating phase transformations in nanostructures requires precise time-integration and long relaxation periods. Users should plan for significant runtime durations.

| Simulation Type | Particle Size | Cores | Approx. Runtime |
| --- | --- | --- | --- |
| **SME (Thermal Cycle)** | 10nm | 8 | ~ 40+ Hours |
| **SME (Thermal Cycle)** | 20nm | 8 | ~ 40+ Hours |
| **SE (Mechanical Load)** | 15nm | 16 | ~ 40+ Hours |

> **Note:** These estimates assume a standard thermal cycling rate (e.g., 20K/ns). Slower heating/cooling rates intended for higher accuracy will linearly increase these times. -->

---

## 🤝 Contributing

Contributions to improve the generator logic or optimize the LAMMPS scripts are welcome.

1. Fork the repository.
2. Create a feature branch (`git checkout -b feature/NewShape`).
3. Commit your changes.
4. Open a Pull Request.

---

## 📚 References & Potentials

The reliability of these simulations depends entirely on the accuracy of the interatomic potential used.

The scripts in this repository are configured to use the **Modified Embedded Atom Method (MEAM)**. Users must ensure that the `NiTi.meam` and `library.meam` files located in the `potentials/` directory correspond to validated literature.

**Recommended Potential:**

> **Ko, W. S., Grabowski, B., & Neugebauer, J. (2015).** Development of a modified embedded-atom method potential for Ni-Ti systems. *Physical Review B*, 92(13), 134107.

*Please ensure you have the appropriate license or rights to use specific potential files if they are not open source.*

---

## 📧 Contact

For questions regarding the workflow implementation or simulation parameters, please open an issue in this repository.

