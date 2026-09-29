# CardioReg assay-identification submission materials

This package contains two manuscript tracks built from the same corrected scientific core.

## AJP-Heart Perspective

`AJP_Perspective/`

- `AJP_Heart_Perspective.pdf` - concise literature-based Perspective to be submitted to the American Journal of Physiology-Heart and Circulatory Physiology.
- `AJP_Heart_Perspective.tex` - source.
- `Cover_Letter.pdf` / `.tex` - submission cover letter.
- `Presubmission_Inquiry.pdf` / `.tex` / `.txt` - optional fit inquiry to the current AJP-Heart editorial office.
- `Submission_Checklist.md` - account-specific items to verify before sending.

The AJP version contains no original unpublished experimental or simulation results. It has no abstract and no figures, consistent with current APS Perspective guidance.

## Full-length manuscript

`Full_Manuscript/`

- `Full_Length_Manuscript.pdf` - detailed framework manuscript.
- `Full_Length_Manuscript.tex` - source.

The full version restores material that is scientifically valid but too detailed for the AJP Perspective: descendant-denominator derivations, transient-versus-cumulative observation operators, shared- versus arm-specific calibration, the full published-data sensitivity calculation, underdetermined structural examples, cluster-level simulations, public-calibration vacuity, signature-invariance stress tests, turnover-model context, and a prospective validation design.

## Reproducibility

- `code/cardioreg_bounds.py` - LP and contrast utilities.
- `code/simulation_utils.py` - biological-unit simulation helpers.
- `code/reproduce_analysis.py` - regenerates all machine-readable analysis tables.
- `tests/test_core.py` - mathematical/software property tests.
- `data/` - generated tables used by the full manuscript.
- `validation/` - source verification, public-calibration audit, and prospective validation plan.
- `REFERENCE_AUDIT.md` - corrected bibliographic records.
- `requirements.txt` and `environment.txt` - computational environment.
- `LICENSE.md` - code and review-material licensing terms.

### Reproduce

From the package root:

```bash
python code/reproduce_analysis.py
pytest -q
```

To rebuild either manuscript, run `pdflatex` twice from its manuscript directory using the corresponding `.tex` file.

## Scope

Neither manuscript claims independent biological validation. The AJP document is a conceptual Perspective. The full manuscript contains synthetic design analyses and a same-study published-data worked example to clarify the framework and its limitations.
