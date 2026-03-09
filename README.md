
# DFVB: Dynamic Flux Vector Balancing

This repository contains the MATLAB code, data, and result files used for the paper "Dynamic balance of sparse flux vectors for efficient simulation of culture dynamics and network reduction". The code implements Dynamic Flux Vector Balancing (DFVB) routines, examples, calibration scripts, validation utilities, and test cases.

## Contents
- Calibration/ — scripts and results used for parameter estimation and bootstrap analyses.
- Functions/ — core DFVB/DFBA helper functions and model builders.
- Models/ — metabolic model files used in examples (MAT files).
- Tests/ — example scripts to run DFVB/DFBA on provided models.
- Validation/ — scripts for experimental comparison, RAS computation and gene activity detection.
- Results/ — output folders and figures produced in the paper.
- README.md — this file.

## Requirements
- MATLAB (R2024a or newer recommended)
- COBRA Toolbox V3.0 or above
- Gurobi

## Quick start
1. Add the repository folder to your MATLAB path.
2. Open and run an example from `Tests/`, e.g. `Tests/Test_ecolicore.m`.
3. Use `Calibration/` scripts to reproduce parameter fitting; result CSV files are provided under `Calibration/Results_Bootstraping/`.
4. See `Validation/` for scripts that compare simulated and experimental data.

## Reproducibility notes
- Each `Results/` subfolder contains Python scripts to reproduce the figures presented in the paper.
- If you plan to run optimization-based routines, ensure required solvers (e.g., Gurobi) and the COBRA Toolbox are installed and configured.

## Citation
If you use this code or data in your work, please cite the paper:

    Tapia García, I., Torrealba, C., Luna, R., Pérez-Correa, J. R., & Saa, P. A. (2026). **Dynamic balance of sparse flux vectors for efficient simulation of culture dynamics and network reduction** 

Contact
For questions contact [iptapia@uc.cl](mailto:iptapia@uc.cl)
