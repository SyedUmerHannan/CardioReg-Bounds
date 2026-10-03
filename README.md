# CardioReg-Bounds

Code, generated tables, and manuscript for **"Target-Specific Partial Identification for Cardiomyocyte Regeneration Assays: Bounds, Sensitivity Analyses, and Design Implications"** (S. U. Hannan).

The repository implements linear-programming bounds on a prespecified target over latent cardiomyocyte histories, given interval-valued assay response probabilities and observation intervals, including endpoint-cell prevalence models that account for death and division.

**Status:** structural/methods work with synthetic design analyses and one worked example on published summary data. No biological validation is claimed. This is not a validated adult-heart regeneration estimator.

## Layout

| Path | Contents |
|---|---|
| `code/cardioreg_bounds.py` | LP bounds (baseline-unit, endpoint-cell, interval signatures), treatment-contrast bounds |
| `code/simulation_utils.py` | Biological-unit (cluster) simulation helpers |
| `code/reproduce_analysis.py` | Regenerates every table in `data/` |
| `tests/test_core.py` | 12 mathematical/software tests |
| `data/` | Generated tables used in the manuscript |
| `validation/` | Source verification, public-calibration audit, prospective validation plan |
| `manuscript/manuscript.tex`, `.pdf` | Manuscript source and build |
| `RESULTS_MAP.md` | Maps each number in the manuscript to the file that produces it |
| `REFERENCE_AUDIT.md` | Corrected bibliographic records |

## Reproduce

```bash
pip install -r requirements.txt
python code/reproduce_analysis.py   # regenerates data/*.csv and analysis_summary.json (seed 20260929)
pytest -q                           # 12 tests
```

Re-running changes only floating-point noise in the last digits of a few values; reported numbers are unchanged. Tested with Python 3.12 / numpy 2.4 / scipy 1.17 and the versions in `environment.txt`.

To rebuild the manuscript, run `pdflatex manuscript.tex` twice in `manuscript/`.

## Scope and limitations

- Worked example: transcription of published summary values from Leone, Musa & Engel, *Cardiovasc Res* 2018 (doi:10.1093/cvr/cvy056). No publisher figure is redistributed. See `validation/Source_Verification.md`.
- Cluster simulations use 100 replicates per cell; coverage values carry Monte Carlo error of up to about 0.05.
- Signature matrices in the synthetic analyses are hypothetical. Public calibration evidence was not found for a generic EdU/pH3/Aurora-B/mononucleation panel (`validation/Public_Calibration_Audit.md`).

## Citation

See `CITATION.cff`. Archived release: Zenodo DOI to be added after the first tagged release.

## License

Code and tests: MIT. Manuscript, figures, and data summaries: see `LICENSE.md`.
