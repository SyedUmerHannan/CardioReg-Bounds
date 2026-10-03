# Results map

Every number in the manuscript comes from the files below. Regenerate with `python code/reproduce_analysis.py`.

| Manuscript item | Value(s) | Source file | Script section |
|---|---|---|---|
| Leone contrast, shared rates, published means (Sec. 4, Eq. 16) | 0.135 to 0.180 | `data/leone_contrast_sensitivity.csv` | `reproduce_analysis.py` (Leone block) |
| Leone contrast, arm-specific rates (Eq. 17) | 0.0245 to 0.2905 | `data/leone_contrast_sensitivity.csv` | same |
| Leone simultaneous Welch difference intervals, shared rates (Eq. 18) | -0.0664 to 0.4049 | `data/leone_contrast_sensitivity.csv`, `data/leone_difference_intervals.csv` | same |
| Published inputs (Fig. 1D values) | means and SDs, n = 5 | `data/leone_published_summary.csv` | same |
| Endpoint prevalence example (Sec. 5) | naive 0.26; correct 0.418; net gain 0.10 | `data/endpoint_denominator_and_net_gain.csv` | endpoint block |
| Net gain bounds without / with count constraint | [-1, 0.294]; [0.05, 0.15] | `data/endpoint_denominator_and_net_gain.csv` | endpoint block |
| Structural bounds, seven histories / four assays (Table 2) | 0.0804 to 0.0940, etc. | `data/structural_target_bounds.csv` | structural block |
| Vacuous bounds with all-[0,1] signatures | 0 to 1 | `data/structural_target_bounds.csv` | structural block |
| Cluster simulation coverage and widths (Sec. 6) | e.g. misspecified coverage 0.62 at n = 30 | `data/cluster_simulation.csv` | simulation block |
| Treatment-dependent signature stress test (Sec. 7) | false-effect fraction 0; widths 0.2 to 0.4 | `data/signature_invariance_null_sensitivity.csv` | invariance block |
| Appendix A propositions | proofs | `manuscript/manuscript.tex`; numerical checks in `tests/test_core.py` | n/a |
