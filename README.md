# ALP-ttbar Toy MC v2.1 — Post-Registration Results: A Likelihood Material Tag in Place of the Hard Material Veto for Displaced Light ALP → e⁺e⁻

**Author:** Won Shik Paik (Independent Researcher, Auckland, New Zealand)
**Date:** 2026-11-12
**Preregistration:** https://doi.org/10.5281/zenodo.24063091
**GitHub (preregistration):** https://github.com/wspaik5/ALP_ttbar_Toy_MC_v2.1_Preregistration_Material_Veto
**Builds on:** v2.0 series overview (https://doi.org/10.5281/zenodo.24056097); v1.7 (https://doi.org/10.5281/zenodo.24001093); v1.8 (https://doi.org/10.5281/zenodo.24034481)
**Motivated by:** arXiv:2610.38932 (Alimena et al.)
**Zenodo (this results record, PDF + zip):** https://doi.org/10.5281/zenodo.24136462
**GitHub (this results record):** https://github.com/wspaik5/ALP_ttbar_Toy_MC_v2.1_Results_Material_Veto

---

## Summary

This repository contains the post-registration results of v2.1, which replaces the hard material veto and opening-angle exemption of v1.1–v1.9 by a single continuous likelihood material tag λ. The tag combines the radial vertex position relative to known material bands and the three-dimensional opening-angle statistic Δχ²_3D in a likelihood ratio. No new simulation or fitting is performed; the tag re-selects the stored, published fits of v1.7 and v1.8.

**Sanity gate:** ALL PASSED (12 published degradations reproduced to machine precision; ∫p_conv = 1 to < 10⁻⁶; λ finite for all base pairs).

**Validations:** T1 (no-material limit), T2 (monotonicity), T3 (resolution dependence) — ALL PASSED.

**Main results (limit degradation from conversions):**

| Config | Mass | Hard veto | Tag | Ratio | Status |
|--------|------|-----------|-----|-------|--------|
| C1 (f=1, pixel z) | 17 MeV | 127.25× | 139.42× | 1.096 | Tag worse |
| C1 | 50 MeV | 63.88× | 65.99× | 1.033 | Tag worse |
| C1 | 100 MeV | 24.24× | 24.52× | 1.011 | ~same |
| C1 | 200 MeV | 6.09× | 6.13× | 1.007 | ~same |
| C2 (f=1, z every layer) | 17 MeV | 117.45× | 111.34× | 0.948 | **Tag better** |
| C2 | 50 MeV | 38.54× | 38.99× | 1.012 | ~same |
| C2 | 100 MeV | 13.14× | 12.80× | 0.974 | **Tag better** |
| C2 | 200 MeV | 3.93× | 3.84× | 0.979 | **Tag better** |
| C3 (f=0.25, z every layer) | 17 MeV | 26.74× | 25.63× | 0.958 | **Tag better** |
| C3 | 50 MeV | 8.83× | 8.73× | 0.989 | ~same |
| C3 | 100 MeV | 3.28× | 2.71× | 0.826 | **Tag better (−17%)** |
| C3 | 200 MeV | 1.09× | 1.09× | 0.997 | ~same |

**Predictions:** 5/7 confirmed (P1, P2, P4, P6, P7), 2/7 refuted (P3, P5).

**Deviations from preregistration:** None.

## Files

- `README.md` — this file
- `v2.1_post_registration_report_2026-11-12.pdf` — full report (12 pages, 3 figures, 2 tables)
- `run_v21.py` — computation script (sanity gate, validations, tag scan, predictions)
- `plot_v21.py` — figure generator (fig15a, fig15b, fig15c)
- `make_report_pdf_v21.py` — report PDF builder
- `fig15a_deg_vs_mass.png` — degradation vs mass (hard veto vs tag)
- `fig15b_lambda_distributions.png` — λ distributions (Config 3, 17 MeV)
- `fig15c_ratio.png` — tag/hard-veto ratio vs mass
- `v21_results.npy` — computed results (numpy dict)

## Checksums

**Report PDF** (`v2.1_post_registration_report_2026-11-12.pdf`):
- SHA-256: `657b7f98258661727400fcdf47d6fdb5887c6c96055713ed2ee985cbdb78d7c2`
- MD5: `dc92692f621540e77c06bd8da5d04404`

## Preregistration integrity

The preregistered README was frozen on 2026-11-09 with:
- SHA-256: `1d625f05f8ee36cf229785aa9cbc747a5f4f031000c27ad595d419e390e35448`
- Zenodo DOI: 10.5281/zenodo.24063091

The tag definition (`mtag_v21.py`), threshold grid, mass windows, counting rule, subsampling and scan logic are exactly as preregistered. The stored fit arrays are loaded unmodified from the published v1.7 and v1.8 results. No post-hoc tuning or additional scans were performed.

## How to reproduce

```bash
# Requires: v1.7 stored fits (v17/*.npz, v17_results.npy),
#           v1.8 stored fits (v18/*.npz, v18_results.npy),
#           v1.9 code (v19/bkg.py, v19/mag.py, v19/run_v13.py, v19/run_v15.py, v19/run_v16.py, v19/run_v17.py),
#           v2.1 preregistration (v21_prereg/mtag_v21.py)
python3 run_v21.py          # runs computation, saves v21_results.npy
python3 plot_v21.py         # generates fig15a, fig15b, fig15c
python3 make_report_pdf_v21.py  # builds report PDF
```

## References

- [E1] CMS Collaboration, JINST 13 (2018) P10034, arXiv:1807.03289.
- [E2] ATLAS Collaboration, JINST 11 (2016) P11020, arXiv:1609.04305.
- [E3] S. F. Koch et al., arXiv:2412.04686.
