Process-Triality Cosmological Framework – Empirical Validation (2025)
==================================================================

This package contains all configuration files, code, and output data
for reproducing the empirical ProTri cosmological run (2025).

Directory structure:
--------------------
Paper/        – Main manuscript (PDF)
Configs/      – CLASS & Cobaya configuration files (YAML + background_protri.c)
Notebook/     – Full execution notebook (Colab-compatible)
Results/      – Posterior chains

How to use:
-----------
1. Open 'Protri_Cobaya_Class_H0onlyFull.ipynb' in Google Colab or Jupyter.
2. Ensure Python 3.10, Cobaya 3.3, and CLASS v2.10 are installed.
3. Run all cells to reproduce the ProTri cosmological results.

Notes:
------
- Chain files are located in Results/chains/
- YAML configurations include active low-ℓ Planck, BAO (6dF), and SN datasets.
- The run converges with R−1 < 0.05 (3000 samples).

Author:
-------
Maximilian Rupp
Independent Researcher, 2025
office@protri.org
