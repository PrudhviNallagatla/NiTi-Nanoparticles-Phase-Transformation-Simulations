# Transformation Exclusion Zone and Stress-Induced Superelastic Response in Free-Standing NiTi Nanoparticles: A Simulation Study

![LAMMPS](https://img.shields.io/badge/LAMMPS-stable__2Aug2023__update3-blue)
![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)
[![DOI](https://img.shields.io/badge/DOI-10.1088/1402--4896/aeb753-blue)](https://doi.org/10.1088/1402-4896/aeb753)

This repository provides LAMMPS input scripts and Python analysis tools developed for the study:

> **Transformation Exclusion Zone and Stress-Induced Superelastic Response in Free-Standing NiTi nanoparticles: A Simulation Approach**  
> Prudhvi Nallagatla, Pradhama Sai Eswar Kumar, Tubin Manuel, Ranjeet Kumar Sahu, Hemantha Kumar, and Debashisha Jena  
> *Physica Scripta*, 2026. [https://doi.org/10.1088/1402-4896/aeb753](https://doi.org/10.1088/1402-4896/aeb753)

## Overview

The provided molecular dynamics (MD) workflows enable researchers to track nanoscale solid-state phase transformations (austenite $\leftrightarrow$ martensite) in equiatomic NiTi nanoparticles. Key capabilities include:
- Modeling size-dependent core–shell transformation dynamics and identifying surface transformation exclusion zones.
- Simulating stress-induced superelasticity (SE) and thermally driven Shape Memory Effects (SME).
- Automated Python routines for atomic trajectory parsing, B2/B19′ phase fractions, and hysteresis characterization.

## Citation

If you use these scripts or data in your research, please cite:

```bibtex
@article{Nallagatla2026NiTi,
  author    = {Nallagatla, Prudhvi and Kumar, Pradhama Sai Eswar and Manuel, Tubin and Sahu, Ranjeet Kumar and Kumar, Hemantha and Jena, Debashisha},
  title     = {Transformation Exclusion Zone and Stress-Induced Superelastic Response in Free-Standing NiTi nanoparticles: A Simulation Approach},
  journal   = {Physica Scripta},
  year      = {2026},
  doi       = {10.1088/1402-4896/aeb753},
  url       = {[https://doi.org/10.1088/1402-4896/aeb753](https://doi.org/10.1088/1402-4896/aeb753)}
}
```

## Theoretical Background

NiTi alloys owe their unique thermomechanical properties to a solid-state phase change between a high-temperature Austenite phase and a low-temperature Martensite phase. 

The Shape Memory Effect occurs when the nanoparticle is mechanically deformed in its Martensite phase. Upon subsequent heating, the material transitions back to Austenite and recovers its original, undeformed shape. 

Superelasticity occurs when mechanical stress is applied while the material is already in the Austenite phase. The applied stress forces a transition into Martensite, allowing the particle to accommodate massive deformation. Upon removal of the load, the material immediately reverts to Austenite and recovers its original shape without thermal activation. In these simulations, this reversible transformation yields a characteristic flag-shaped hysteresis loop.

## Computational Methodology

The simulation workflow requires sequential execution. Before conducting mechanical deformation tests, a thermal baseline must be established. This initial thermal cycling step equilibrates the nanoparticle and triggers the baseline phase transformation.

Following the thermal cycling, use the provided Python scripts to validate the phase fractions. The nanoparticle must be 100% Martensite at the target low temperature and 100% Austenite at the target high temperature. Once these thermal baselines are confirmed, the mechanical tests for SME and SE can proceed.

## Installation & Dependencies

To ensure a reproducible environment and avoid local dependency conflicts, we recommend running this repository either as a Docker container or as a VS Code Devcontainer. 

The included setup script (`all_dependencies.sh`) handles the installation of all required components: system packages, Python libraries (NumPy, SciPy, Matplotlib, pandas, ase, ovito), and LAMMPS `stable_2Aug2023_update3` configured with the necessary packages (MEAM, MANYBODY, PTM, VORONOI, and MPI).

To set up the environment manually within a container, run:

```bash
chmod +x all_dependencies.sh
./all_dependencies.sh
```

## Execution Protocol

Execute the simulations in the following sequence:

1. **Establish the Thermal Baseline:**
   ```bash
   lmp -in src/lmps/thermal_cycle.lmp
   ```

2. **Validate Phase Fractions:**
   Use the Python scripts in `scripts/analysers/` to check the composition. The nanoparticle should be 100% Martensite at the low-temperature limit and 100% Austenite at the high-temperature limit.

3. **Execute Mechanical Tests:**
   Once the thermal baseline is validated, you can run the mechanical tests concurrently.
   
   For the Shape Memory Effect:
   ```bash
   lmp -in src/lmps/sme.lmp
   ```
   
   For Superelasticity:
   ```bash
   lmp -in src/lmps/se.lmp
   ```

## Repository Architecture

Here is how the repository is organized:

* `src/lmps/`: LAMMPS input scripts for energy minimization, thermal cycling, SME, and SE simulations.
* `scripts/`: Python post-processing tools for hysteresis loop analysis, phase fraction calculations, and figure generation.
* `potentials/`: Interatomic potential files (MEAM potential for Ni-Ti systems).
* `data/`: Initial crystallographic structures and simulation output trajectories.
* `paper/` and `docs/`: LaTeX manuscript files and supplementary documentation.

## Data Availability

Once the associated paper is published, we will upload the full simulation dataset and trajectories to the [NOMAD](https://nomad-lab.eu/) repository and add a permanent DOI to this README.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for complete details.

## References

* Ko, W. S., Grabowski, B., & Neugebauer, J. (2015). Development of a modified embedded-atom method potential for Ni-Ti systems. *Physical Review B*, 92(13), 134107.
* Waitz, T., Antretter, T., Keckes, J., & Karnthaler, H. P. (2004). Size-dependent martensitic transformation in NiTi nanoparticles. *Acta Materialia*, 52(1), 137-147.
