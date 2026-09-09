GPNMB UKB-PPP OLINK -> ieu-b-7 LD-CLUMPED MR RESULTS SUMMARY
============================================================

Purpose
-------
This file summarises the final Olink/UKB-PPP marginal-data Mendelian
randomisation results stored in this directory.

The results below are read from:
- 04_all_ld_clumped_mr_results.csv
- 05_multiplicative_random_effects_ivw_results.csv

The saved LD-clumped instrument files are the final instrument sets used
for these calculations. No new OpenGWAS API request is required to
reproduce the random-effects calculations from the saved clumps.

Input/QC summary
----------------
Regional harmonised exposure-outcome records: 4,338
Records retained after QC + genome-wide significance filtering: 790

QC included:
- INFO >= 0.8
- allele frequency 0.01-0.99
- F-statistic >= 10
- exposure p < 5e-8
- outcome N >= 450,000

LD-clumping specifications
--------------------------
- strict:      r2 < 0.001 within 10 Mb
- standard:    r2 < 0.01 within 10 Mb
- sensitivity: r2 < 0.01 within 1 Mb
- liberal:     r2 < 0.05 within 1 Mb

Retained instruments
--------------------
- strict: 5
- standard: 12
- sensitivity: 12
- liberal: 25


FIXED-EFFECT IVW RESULTS
========================

Setting                          n       Beta         SE           p        OR                  95% CI          Q  Q df         Q p
-----------------------------------------------------------------------------------------------------------------------------------
strict_r2_0p001_10Mb             5   0.075957   0.039869   0.0567574    1.0789           0.9978-1.1666    28.2085     4   1.132e-05
standard_r2_0p01_10Mb           12   0.089754   0.037437   0.0165087    1.0939           1.0165-1.1772    42.2158    11   1.486e-05
sensitivity_r2_0p01_1Mb         12   0.089754   0.037437   0.0165087    1.0939           1.0165-1.1772    42.2158    11   1.486e-05
liberal_r2_0p05_1Mb             25   0.105485   0.032850  0.00132213    1.1112           1.0420-1.1852    61.6171    24   3.763e-05


MULTIPLICATIVE RANDOM-EFFECTS IVW RESULTS
=========================================

Setting                          n       Beta      RE SE        RE p     RE OR               RE 95% CI       phi
----------------------------------------------------------------------------------------------------------------
strict_r2_0p001_10Mb             5   0.075957   0.105875    0.473113    1.0789           0.8767-1.3277    7.0521
standard_r2_0p01_10Mb           12   0.089754   0.073340    0.221026    1.0939           0.9474-1.2630    3.8378
sensitivity_r2_0p01_1Mb         12   0.089754   0.073340    0.221026    1.0939           0.9474-1.2630    3.8378
liberal_r2_0p05_1Mb             25   0.105485   0.052636   0.0450624    1.1112           1.0023-1.2320    2.5674


Interpretation
--------------
The strict and standard clumping specifications are the principal
multi-instrument analyses.

Fixed-effect IVW:
- Strict clumping is borderline positive:
  beta = 0.075957,
  OR = 1.0789,
  p = 0.056757.

- Standard clumping is nominally positive:
  beta = 0.089754,
  OR = 1.0939,
  p = 0.016509.

Heterogeneity:
- Cochran's Q is highly significant under both principal settings,
  indicating substantial between-instrument heterogeneity.

Multiplicative random-effects IVW:
- Strict:
  beta = 0.075957,
  SE = 0.105875,
  OR = 1.0789,
  95% CI 0.8767-
  1.3277,
  p = 0.473113.

- Standard:
  beta = 0.089754,
  SE = 0.073340,
  OR = 1.0939,
  95% CI 0.9474-
  1.2630,
  p = 0.221026.

The multiplicative random-effects model retains the fixed-effect IVW
point estimate but inflates its standard error according to residual
heterogeneity. The strict and standard random-effects results are
non-significant and therefore temper interpretation of the corresponding
fixed-effect estimates.

The liberal setting remains sensitivity evidence only. Its nominal
random-effects significance should not replace the principal strict and
standard analyses.

Reproducibility
---------------
Fresh LD clumping through run_gpnmb_ld_clumped_mr.R requires valid
OpenGWAS authentication because it requests European-reference LD
clumping.

The supplied script:

01_scripts/run_gpnmb_random_effects_from_saved_clumps.R

recalculates the fixed-effect IVW estimates, Cochran's Q and
multiplicative random-effects IVW estimates directly from the saved final
clumped instrument files in this directory and does not require a new
OpenGWAS API request.

Key output files
----------------
- 01_filtered_genomewide_significant_harmonised_variants.csv
- 02_clumped_instruments_strict_r2_0p001_10Mb.csv
- 02_clumped_instruments_standard_r2_0p01_10Mb.csv
- 02_clumped_instruments_sensitivity_r2_0p01_1Mb.csv
- 02_clumped_instruments_liberal_r2_0p05_1Mb.csv
- 03_mr_result_*.csv
- 04_per_variant_ratios_*.csv
- 04_all_ld_clumped_mr_results.csv
- 05_multiplicative_random_effects_ivw_results.csv

Interpretation rule
-------------------
Do not treat all regional variants as statistically independent
instruments. Strict and standard LD-clumped analyses are the principal
multi-instrument Olink results; liberal clumping is sensitivity only.
