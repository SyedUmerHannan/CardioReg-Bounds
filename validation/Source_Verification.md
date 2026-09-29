# Source verification: Leone, Musa, and Engel 2018 worked example

Primary source:

Leone M, Musa G, Engel FB. *Cardiomyocyte binucleation is associated with aberrant mitotic microtubule distribution, mislocalization of RhoA and IQGAP3, as well as defective actomyosin ring anchorage and cleavage furrow ingression.* Cardiovascular Research. 2018;114(8):1115-1131. DOI: 10.1093/cvr/cvy056. PMID: 29522098.

Publisher locations used:

- **Methods, section 2.5, "Statistical analysis"**: states that data are expressed as mean +/- standard deviation (SD) from independent experiments.
- **Figure 1D and its caption**: reports two-sided, one-sided, and no-furrow ingression for P3 cardiomyocytes stimulated with FBS or FGF1/p38 inhibition; the caption states `n > 162 cells from five independent experiments` and `Data are mean +/- SD`.
- **Figure 1E and Results section 3.1**: reports that >90% of P3 cardiomyocytes with two-sided ingression divided into two mononucleated daughter cells, whereas >85% of cells with one-sided or no furrow ingression failed cytokinesis.
- **Results section 3.1**: gives the FBS and FGF1/p38i morphology means used in the package: FBS 25 +/- 13%, 33 +/- 14%, 42 +/- 9%; FGF1/p38i 43 +/- 10%, 26 +/- 13%, 31 +/- 5%.

The package does not redistribute the publisher's figure. `data/leone_published_summary.csv` is a transcription of published summary values.

Important limitation: the prose outcome thresholds (>90% and >85%) do not provide category-specific denominators sufficient to treat them as statistical confidence limits. The manuscript therefore uses them only as sensitivity ranges.
