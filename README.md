# GPNMB PARKINSON'S DISEASE DISSERTATION — SUPPLEMENTARY PACKAGE

## Purpose
This supplementary package accompanies the dissertation:

"Prioritising GPNMB as a Parkinson's disease-relevant protein using
Mendelian randomisation and colocalisation"

It contains the analysis scripts, derived results, selected supporting
inputs and reproducibility records used to support the reported GPNMB
Mendelian randomisation, colocalisation, SuSiE and cross-tissue analyses.

## PACKAGE STRUCTURE

### 01_scripts/
Analysis and reproducibility scripts.

Key files include:
- GPNMB_MR_full_pipeline_with_olink_rsid_mapping_v2.R
  Open Targets Olink variant formatting, rsID mapping and harmonisation.

- GPNMB_coloc_peer_review_pipeline.R
  Main regional coloc.abf workflow.

- GPNMB_MAF_check_and_coloc_susie.R
  Regional coverage/MAF checks and LD-aware SuSiE sensitivity analysis.

- run_gpnmb_clumping_lead_mr.R
  Lead-variant Olink MR analysis.

- run_gpnmb_ld_clumped_mr.R
  QC, OpenGWAS European-reference LD clumping, fixed-effect IVW and
  multiplicative random-effects IVW workflow. A fresh clumping run
  requires valid OpenGWAS authentication.

- run_gpnmb_random_effects_from_saved_clumps.R
  Recalculates fixed-effect IVW, Cochran's Q and multiplicative
  random-effects IVW directly from the saved final clumped instrument
  sets. This script does not require a new OpenGWAS API request.

- rs5850_independent_pd_lookup_v2.R
  Supporting rs5850 Parkinson's disease outcome lookup, including proxy
  search where required.

- generate_figures.R
  Reproduces Figures 2-5 from the final saved analytical outputs using
  base R graphics. Figures 1a and 1b were assembled manually in
  diagrams.net (draw.io).


### 02_results/
Derived analytical results.

02_results/current_results/
  SomaScan MR outputs, regional colocalisation inputs and results,
  SuSiE outputs, QC/audit files, Open Targets mapping/harmonisation
  records and transcriptomic supporting summaries.

Important examples:
- somascan/somascan_gpnmb_pd_odds_ratios.csv
- gpnmb_final_coloc_results_peer_review.csv
- gpnmb_coloc_prior_sensitivity.csv
- gpnmb_coloc_susie_results.csv
- qc/gpnmb_coloc_coverage_audit.csv
- transcriptomics/gpnmb_full_eqtl_coloc_summary.csv

02_results/latest_olink/
  Final UKB-PPP Olink marginal-data MR pathway.

Key files:
- gpnmb_ukbppp_ieu_b_7_harmonised_region_overlap.csv
  4,338 harmonised regional exposure-outcome records before final QC.

- gpnmb_corrected_mr_outputs/01_lead_variant_wald_mr.csv
  Lead-variant rs75801644 Wald-ratio result.

- gpnmb_ld_clumped_mr_outputs/
  Contains the 790 QC-filtered genome-wide significant records,
  saved strict/standard/sensitivity/liberal clumped instrument sets,
  per-variant ratios, fixed-effect IVW results and multiplicative
  random-effects IVW results.

Principal files in this folder:
- 01_filtered_genomewide_significant_harmonised_variants.csv
- 04_all_ld_clumped_mr_results.csv
- 05_multiplicative_random_effects_ivw_results.csv
- 00_READ_ME_RESULTS_SUMMARY.txt

The strict and standard clumping specifications are the principal
multi-instrument Olink analyses. The liberal specification is retained
as sensitivity evidence only.


### 03_supporting_inputs/
Selected source, annotation and molecular-QTL files used to support or
audit the analyses.

Contents include:
- prot-c-5080_131_3.vcf.gz
  Original SomaScan GPNMB exposure VCF used for the regional coverage
  audit and colocalisation source-level checks.

- OpenTargets/ENSG00000136235-qtl-credible-sets-target.tsv
  GPNMB Open Targets molecular-QTL credible-set export.

- GTEx/
  rs5850 eQTL/sQTL and related GTEx exports.

- expression_atlas/
  Cross-tissue expression-support files.

- ensembl_rs5850/
  Variant annotation, overlapping-feature and LD-proxy records.

- gwas_catalog_rs5850/
  GWAS Catalog-style rs5850 association/study/trait exports.

- independent_pd_outcomes/
  Metadata used for independent Parkinson's disease outcome checks.

- olink_ukbppp_regional/
  Selected regional UKB-PPP exposure and Parkinson's disease outcome
  extracts plus the Olink protein metadata map:
  * gpnmb_ukbppp_chr7_22p8_23p8Mb_exposure.tsv
  * ieu_b_7_chr7_22p8_23p8Mb_PD_outcome.vcf
  * olink_protein_map_3k_v1.tsv


### 04_supporting_evidence/
Supporting cross-layer biological and molecular evidence.

Files include:
- EMS83968-supplement-Supplementary_tables.xlsx
  Published supplementary material containing the GPNMB.5080.131.3
  SomaScan assay and associated pQTL/eQTL/disease-overlap information.

- gpnmb_brain_eqtl_top_snp_summary.csv
  Brain eQTL summary for key GPNMB variants.

- gpnmb_eqtl_pqtl_key_snp_overlap.csv
  Cross-layer eQTL/pQTL overlap summary.


### 05_reproducibility/
Software-session, figure-generation and source-provenance records.

Files include:
- MR_sessionInfo.txt
  R/package environment used for Mendelian randomisation analyses.

- GPNMB_coloc_sessionInfo.txt
  R/coloc/susieR environment used for colocalisation analyses.

- figure_generation_notes.txt
  Documents how Figures 1a-5 were prepared, identifies the source files
  and values used for Figures 2-5, records their reproduction using base
  R graphics, and documents that Figures 1a-1b were assembled manually
  in diagrams.net (draw.io).

- source_provenance/
  Source-file and regional-processing provenance:
  * 00_file_info.txt
  * gpnmb_processed_summary.txt


## ANALYTICAL FLOW

### SomaScan MR
SomaScan GPNMB exposure -> Parkinson's disease outcome harmonisation ->
single-instrument rs5850 Wald ratio.

Because only one valid instrument was available, the SomaScan result is
supportive single-SNP evidence and cannot provide multi-instrument
heterogeneity or horizontal-pleiotropy diagnostics.


### UKB-PPP Olink MR
UKB-PPP chromosome 7 GPNMB marginal summary data + ieu-b-7 ->
4,338 harmonised regional records -> QC and genome-wide significance ->
790 QC-filtered records -> rsID deduplication before European-reference LD
clumping -> strict, standard,
sensitivity and liberal instrument sets -> fixed-effect IVW.

Cochran's Q showed substantial heterogeneity. Multiplicative
random-effects IVW was therefore calculated as a sensitivity analysis
from the saved final clumped instrument sets. The strict and standard
random-effects results are non-significant and are used to temper the
fixed-effect interpretation.


### Regional colocalisation
Original SomaScan regional source: 165 variants in the analysed interval.

After overlap/harmonisation:
- coloc.abf primary analysis: 150 SNPs
- SuSiE sensitivity analysis: 146 SNPs

The primary coloc.abf result supports a shared signal under the main
prior, while prior-sensitivity and SuSiE results show that interpretation
depends on prior assumptions and unresolved multi-signal locus structure.


## REPRODUCIBILITY NOTES

- The saved LD-clumped Olink instrument files are the final instrument
  sets used for the reported analyses.

- A new OpenGWAS LD-clumping request requires valid OpenGWAS
  authentication. The supplied
  run_gpnmb_random_effects_from_saved_clumps.R script reproduces the
  fixed-effect IVW, Cochran's Q and multiplicative random-effects IVW
  calculations directly from the saved clumped sets without making a
  new API request.

- The full very large third-party source datasets are not redistributed
  in this package. Selected regional extracts, derived harmonised files,
  source identifiers and provenance records are supplied instead.

- Figures 2-5 can be reproduced from the retained final outputs using
  01_scripts/generate_figures.R. Figure-generation provenance is
  documented in 05_reproducibility/figure_generation_notes.txt.

- Open Targets, OpenGWAS, UKB-PPP/Synapse, GTEx, Expression Atlas,
  Ensembl and other third-party resources remain subject to the terms
  and citation requirements of their original providers.

- The dissertation should be used for the formal scientific
  interpretation of these outputs. The supplementary package is intended
  to document and support reproducibility of the reported analyses.


## MAIN RESULTS TO CROSS-CHECK

Olink:
- Regional harmonised records: 4,338
- QC-filtered genome-wide significant records: 790
- Strict clumping: 5 instruments
- Standard clumping: 12 instruments
- Sensitivity clumping: 12 instruments
- Liberal clumping: 25 instruments

Principal fixed-effect IVW:
- Strict: OR approximately 1.08; p approximately 0.057
- Standard: OR approximately 1.09; p approximately 0.017

Principal multiplicative random-effects IVW:
- Strict: OR approximately 1.08; p approximately 0.473
- Standard: OR approximately 1.09; p approximately 0.221

Colocalisation:
- Primary coloc.abf input: 150 SNPs
- SuSiE input: 146 SNPs
- Main-prior PP.H4: approximately 0.944

These values should agree with the dissertation tables and the
machine-readable result files included in this package.
