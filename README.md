# Transformation Exclusion Zone and Stress-Induced Superelastic Response in Free-Standing NiTi Nanoparticles: A Simulation Study

![LAMMPS](https://img.shields.io/badge/LAMMPS-stable__2Aug2023__update3-blue)
![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![DOI](https://img.shields.io/badge/DOI-Pending-lightgrey)

## Overview

This repository provides a comprehensive computational framework for simulating the thermomechanical behaviors of Nickel-Titanium (NiTi) nanoparticles using Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS). The developed input scripts and integrated Python analysis tools enable the precise modeling of complex phase transformations at the nanoscale. By leveraging this workspace, researchers can systematically investigate both the Shape Memory Effect (SME) and Superelasticity (SE) phenomena, facilitating a deeper understanding of martensitic transformations in advanced shape memory alloys.

## Theoretical Background

Nickel-Titanium alloys exhibit remarkable thermomechanical properties driven by solid-state phase transformations between the high-temperature Austenite phase and the low-temperature Martensite phase. The Shape Memory Effect manifests when the nanoparticle undergoes mechanical deformation while in the Martensite phase. Upon subsequent heating, the material transitions back to the Austenite phase, entirely recovering its original, undeformed geometry. This temperature-induced structural recovery forms the basis of numerous engineering and biomedical applications.

In contrast, Superelasticity occurs when mechanical stress is applied to the nanoparticle while it is in the high-temperature Austenite phase. The applied stress induces a transition to the Martensite phase, allowing the material to accommodate massive deformations. Upon the removal of the mechanical load, the material spontaneously reverts to the Austenite phase and recovers its original shape without requiring thermal activation. This reversible, stress-induced transformation yields a characteristic flag-shaped hysteresis loop during the loading and unloading cycles.

## Computational Methodology

The simulation workflow requires a strict sequential execution to ensure thermodynamic stability and accurate phase representations. Initially, a thermal baseline must be established to equilibrate the nanoparticle and facilitate the baseline phase transformation. Following this thermal cycling, the phase fractions must be rigorously validated. It is critical to confirm that the nanoparticle consists entirely of the Martensite phase at the target low temperature and the Austenite phase at the target high temperature. Only after verifying these thermal baselines can the subsequent mechanical tests for the Shape Memory Effect and Superelasticity be executed accurately.

## Installation & Dependencies

To ensure a fully reproducible computational environment and avoid local dependency conflicts, it is highly recommended to run this repository either as a standalone Docker container or as a Devcontainer within Visual Studio Code. 

A comprehensive bash script (`all_dependencies.sh`) is provided for environment configuration. When building the container, this script automatically handles the installation of all required dependencies, including system packages, Python libraries (such as NumPy, SciPy, Matplotlib, pandas, ase, and ovito), and the specific version of LAMMPS (`stable_2Aug2023_update3`) configured with necessary packages (MEAM, MANYBODY, PTM, VORONOI, and MPI).

If setting up the environment manually within a container, run:

```bash
chmod +x all_dependencies.sh
./all_dependencies.sh
```

## Execution Protocol

The simulations must be executed sequentially, establishing the thermal baseline before proceeding to the mechanical deformation tests.

1. **Establish the Thermal Baseline:**
   ```bash
   lmp -in src/lmps/thermal_cycle.lmp
   ```

2. **Validate Phase Fractions:**
   Utilize the Python analysis scripts located in the `scripts/analysers/` directory. Ensure the nanoparticle exhibits a 100% Martensite composition at the low-temperature limit and a 100% Austenite composition at the high-temperature limit.

3. **Execute Mechanical Tests:**
   Once the thermal baseline is validated, the mechanical tests can be run concurrently.
   
   For the Shape Memory Effect:
   ```bash
   lmp -in src/lmps/sme.lmp
   ```
   
   For Superelasticity:
   ```bash
   lmp -in src/lmps/se.lmp
   ```

## Repository Architecture

The workspace is organized to separate simulation inputs, analytical tools, and documentation:

* `src/lmps/`: LAMMPS input scripts encompassing energy minimization, thermal cycling, Shape Memory Effect, and Superelasticity simulations.
* `scripts/`: Python post-processing utilities designed for hysteresis loop analysis, phase fraction calculation, and publication-quality figure generation.
* `potentials/`: Interatomic potential files, specifically the modified embedded-atom method (MEAM) potential for Ni-Ti systems.
* `data/`: Initial crystallographic structures and resulting simulation output trajectories.
* `paper/` and `docs/`: LaTeX manuscript source files and supplementary documentation.

## Data Availability

Upon publication, the full simulation dataset and trajectories from this repository will be archived and openly accessible on the [NOMAD (Novel Materials Discovery)](https://nomad-lab.eu/) repository, at which point a permanent DOI will be minted and added to this document.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for complete details.

## References

* Ko, W. S., Grabowski, B., & Neugebauer, J. (2015). Development of a modified embedded-atom method potential for Ni-Ti systems. *Physical Review B*, 92(13), 134107.
* Waitz, T., Antretter, T., Keckes, J., & Karnthaler, H. P. (2004). Size-dependent martensitic transformation in NiTi nanoparticles. *Acta Materialia*, 52(1), 137-147.
