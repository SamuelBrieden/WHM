# WHMcode — figures of the web-halo model paper

Jupyter notebooks that produce the figures of Brieden et al. (2025), *The web-halo model: a halo model for the cosmic web*
(arXiv:2508.10902), from scratch with the modified CAMB of this repository (`../WHM-CAMB`).

| notebook | figures |
|---|---|
| `WHMcode.ipynb` | Fig. 2 `halomassfunctions.png`, Fig. 4 `WHMcode_spectra_comp_kNL1_concnfw_narya.png` (+`_leg_wide`), Fig. 5 `w0waMnu_LPT_kNL1.0_BACCO_baselineaccuratebf_concnfwwarren_100.png`, Fig. 7 `theo_vs_measured_p1h_linlag_all.png` (+`_leg_v2`), Fig. 8 `HMF_cumulative_integrand_comparison.png`, Fig. A1 `test_cutoff_and_jn.png` |
| `whmcode_windowfunctions.ipynb` | Fig. 3 `WHMcode_wk_v3.png` (analytic, self-contained) |
| `WHMcode_LPT_scatter.ipynb` | Fig. 6 `EmulatorComparison.png` (+`_leg_wide`), Fig. A2 `LPT_scatter.png` (plotting only; reads the per-cosmology spectra written by `WHMcode.ipynb`) |

Fig. 1 and the appendix images of the window-function measurements (`LPT-Filament-Halo-density.png`, `8haloes_per_filament*.png`,
`pk_10Mpc_2_mod_correct_8haloes.png`) are N-body visualisations and are not produced here.

## Inputs (`data/`)

- `Pk_Marcos/`: Planck-2018 N-body measurements (linear/non-linear P(k), 1-halo term, halo mass function, u(k,M)) used in Figs. 2, 7, 8.
- `CosmicEmu/`: the 112 Mira-Titan cosmologies (`MiraTitanParams.dat`) and their emulator spectra at z = 0, 0.8, 1.5 (Fig. 6, third row).
- `cosmo/`: the three cosmology samples (`*_w0wamnu_BACCO.txt`, `*_w0wamnu_EE2.txt`, `*_w0wamnu_CosmicEmu.txt`; 100, 100 and 112 models).
  Only physical densities (`om_cdm`, `om_b`, `om_nu`) are stored; the fractional ones are derived in the notebook.

Everything else under `data/` (`Pklin/`, `PkHM/`, `PkHMori/`, `PkWHMhf/`, `PkLPT/`, `PkBacco/`, `PkEE2/`, `knl/`, `fiducial_bacco/`,
`Pk_Marcos/*_warren.txt`) is a cache written by `WHMcode.ipynb` when `WRITE` is enabled and reloaded on later runs.

## Running

`WHMcode.ipynb` runs top to bottom from this directory. The configuration cell defines

- `RECOMPUTE[...]` / `WRITE[...]` per input source (`camb`, `velocileptors`, `bacco`, `ee2`, `miratitan`, `nbody`):
  with everything `False` the cached spectra under `data/` are used; with everything `True` all of them are regenerated.
- `SAMPLE`: the cosmology sample of sections 4-7, `"BACCO"` (paper default; Fig. 5 and the first row of Fig. 6), `"EE2"` or `"CosmicEmu"` (rows 2 and 3 of Fig. 6).

For headless runs the same settings can be given as environment variables, e.g.

```bash
WHM_SAMPLE=BACCO WHM_RECOMPUTE=all WHM_WRITE=all jupyter nbconvert --to notebook --execute WHMcode.ipynb
WHM_SAMPLE=EE2   WHM_RECOMPUTE=all WHM_WRITE=all jupyter nbconvert --to notebook --execute WHMcode.ipynb
WHM_SAMPLE=CosmicEmu WHM_RECOMPUTE=all WHM_WRITE=all jupyter nbconvert --to notebook --execute WHMcode.ipynb
jupyter nbconvert --to notebook --execute WHMcode_LPT_scatter.ipynb
jupyter nbconvert --to notebook --execute whmcode_windowfunctions.ipynb
```

(`WHM_NCOSMO` overrides the number of sample cosmologies.) A full regeneration takes several hours, dominated by the
per-cosmology CAMB and velocileptors runs; the cached mode takes about ten minutes.

## Model variants used by the figures

`WHM-CAMB` exposes the web-halo model through `halofit_version='brieden2025_*'` plus two switches passed to
`NonLinearModel.set_params`: `WHM_thinweb=1` selects the 'thin web' sheet/filament profiles (upper edges of the coloured bands
in Fig. 4), and `WHM_hmf` selects the halo mass function: 0 Sheth & Tormen (baseline, used everywhere unless stated),
1 Warren et al. 2006 (fourth row of Fig. 5 and the purple curves of Fig. 8), 2 Despali et al. 2016 (available, not used by any figure).

## velocileptors grid

The 1-loop LPT spectra are computed with velocileptors' `LPT_RSD` (`k_IR = k_cutoff = sqrt5 k_nl(z)`, `jn=15`). The paper used
velocileptors as of July 2025 (commit 83d82ae), whose FFTLog grid defaults were `N=2000, extrap_min=-5, extrap_max=3`. Commit
f8dae25 (2026-04-23, "Replaced old module", released as velocileptors 3.1, the version pinned here) changed those defaults to
`N=256, extrap_min=-4, extrap_max=1`, which alters the spectra by 1-3 % below k = 1 h/Mpc and by up to 20 % of P_NL beyond the
cutoff. The notebook therefore passes the July-2025 values explicitly (`LPT_GRID` in the configuration cell); with them
velocileptors 3.1 reproduces the paper's spectra to better than 1e-3.

## Environment

`env/environment.yml` is the pinned conda environment (`conda env create -f env/environment.yml`), after which the local
`WHM-CAMB` is built with `pip install -e ../WHM-CAMB --no-build-isolation` and `classy` from the public CLASS 3.3.4
(https://github.com/lesgourg/class_public, tag v3.3.4: `make -j4 class; make libclass.a; CC=gcc pip install . --no-build-isolation`;
the notebook only needs standard classy and adds the HMcode displacement scale itself); `env/import_smoketest.py` checks the imports.
On macOS/arm64 a freshly built `camblib.so` may need `codesign --force --sign - camb/camblib.so` before it can be loaded.
The figures use matplotlib with `usetex`, so a LaTeX installation is required.

## AI Statement

The original code, plotting scripts, and paper were created by me without any LLM assistance. 
However, I used Claude Code to clean the notebooks, document the environment, clarify the user flags in the CAMB fortran source code, and generate this README file. 
Claude also helped me finding out that the discrepancy with respect to the published papers originated from the change in velocileptors specified and fixed above. 
I verified all new additions from Claude and take full responsibility for their validity. 
