# Public calibration audit

The literature review supporting this package did not identify a public cardiomyocyte-regeneration dataset that provides a transportable conditional-response matrix for the same terminal-history partition across a generic EdU/pH3/Aurora-B/mononucleation panel.

Accordingly, the entries in `public_calibration_audit.csv` remain `[0,1]` when direct calibration evidence is unsupported. In the synthetic structural check, replacing the hypothetical signature matrix with all `[0,1]` entries makes the productive-division target exactly `[0,1]`. This is intentional: the framework does not manufacture information when assay-to-history calibration is absent.

The audit is not a fitted signature matrix and should not be used as biological calibration.
