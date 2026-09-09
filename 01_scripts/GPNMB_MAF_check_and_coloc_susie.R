# GPNMB_MAF_check_and_coloc_susie.R
#
# Purpose:
#   1. Check whether OpenGWAS returns eaf/MAF for the GPNMB pQTL and PD regional data.
#   2. If pQTL eaf/MAF is missing, document that explicitly.
#   3. Build an LD matrix for the harmonised coloc SNPs using ieugwasr::ld_matrix().
#   4. Run coloc.susie as an LD-aware multi-signal sensitivity analysis.
#
# Inputs expected:
#   - gpnmb_coloc_merged_harmonised_peer_review.csv
#
# Outputs:
#   - gpnmb_opengwas_raw_maf_eaf_qc.csv
#   - gpnmb_coloc_susie_ld_matrix.csv
#   - gpnmb_coloc_susie_input_variants.csv
#   - gpnmb_coloc_susie_exposure_summary.csv
#   - gpnmb_coloc_susie_outcome_summary.csv
#   - gpnmb_coloc_susie_summary.csv
#   - gpnmb_coloc_susie_results.csv
#
# Caveats:
#   - OpenGWAS LD reference may drop variants absent from the LD panel.
#   - The OpenGWAS LD panel is a reference panel, not your exact study genotype data.
#   - pQTL-specific allele frequencies may remain unavailable from the API.
#   - sdY = 1 assumes the pQTL phenotype is standardised.

library(tidyverse)
library(ieugwasr)
library(coloc)

if (!requireNamespace("susieR", quietly = TRUE)) {
  install.packages("susieR")
}
library(susieR)

# Allow running either from the package root or from 01_scripts.
if (basename(normalizePath(getwd(), winslash = "/", mustWork = FALSE)) == "01_scripts") {
  setwd("..")
}

results_dir <- file.path("02_results", "current_results")
susie_dir <- file.path(results_dir, "susie")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(susie_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 0. Parameters
# ==============================================================================

chr_locus <- "7"
start_pos <- 22814547
end_pos   <- 23814547
region_query <- paste0(chr_locus, ":", start_pos, "-", end_pos)

pqtl_id <- "prot-c-5080_131_3"
pd_id   <- "ieu-b-7"

pqtl_N <- 996

pd_ncase <- 33674
pd_ncontrol <- 449056
pd_N <- pd_ncase + pd_ncontrol
pd_s <- pd_ncase / pd_N

merged_path <- file.path(results_dir, "gpnmb_coloc_merged_harmonised_peer_review.csv")

# ==============================================================================
# 1. Helpers
# ==============================================================================

safe_nrow <- function(x) {
  if (is.null(x)) 0 else nrow(x)
}

pick_col <- function(dat, candidates) {
  hit <- candidates[candidates %in% names(dat)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

summarise_maf_eaf_from_raw <- function(dat, label) {
  if (is.null(dat) || nrow(dat) == 0) {
    return(tibble(
      dataset = label,
      raw_rows = 0,
      columns = NA_character_,
      eaf_column = NA_character_,
      maf_column = NA_character_,
      eaf_nonmissing = 0,
      maf_nonmissing = 0
    ))
  }

  dat <- as_tibble(dat)

  eaf_col <- pick_col(dat, c("eaf", "EAF", "af", "AF", "freq"))
  maf_col <- pick_col(dat, c("maf", "MAF"))

  eaf_nonmissing <- if (!is.na(eaf_col)) sum(!is.na(dat[[eaf_col]])) else 0
  maf_nonmissing <- if (!is.na(maf_col)) sum(!is.na(dat[[maf_col]])) else 0

  tibble(
    dataset = label,
    raw_rows = nrow(dat),
    columns = paste(names(dat), collapse = ", "),
    eaf_column = eaf_col,
    maf_column = maf_col,
    eaf_nonmissing = eaf_nonmissing,
    maf_nonmissing = maf_nonmissing
  )
}

download_raw_region <- function(id, label) {
  message("Downloading raw OpenGWAS region for ", label, ": ", id)

  tryCatch(
    {
      ieugwasr::associations(
        variants = region_query,
        id = id,
        proxies = 0
      )
    },
    error = function(e) {
      message("Download failed for ", label, ": ", conditionMessage(e))
      NULL
    }
  )
}

# ==============================================================================
# 2. Check whether OpenGWAS exposes eaf/MAF for pQTL and PD
# ==============================================================================

pqtl_raw <- download_raw_region(pqtl_id, "pQTL")
pd_raw   <- download_raw_region(pd_id, "PD GWAS")

maf_qc <- bind_rows(
  summarise_maf_eaf_from_raw(pqtl_raw, "GPNMB pQTL OpenGWAS raw"),
  summarise_maf_eaf_from_raw(pd_raw, "PD GWAS OpenGWAS raw")
)

print(maf_qc)
write_csv(maf_qc, file.path(results_dir, "gpnmb_opengwas_raw_maf_eaf_qc.csv"))

if (maf_qc$eaf_nonmissing[maf_qc$dataset == "GPNMB pQTL OpenGWAS raw"] == 0 &&
    maf_qc$maf_nonmissing[maf_qc$dataset == "GPNMB pQTL OpenGWAS raw"] == 0) {
  message("pQTL-specific eaf/MAF was not returned by OpenGWAS. You need original pQTL summary stats to obtain pQTL-specific allele frequencies.")
}

# ==============================================================================
# 3. Load harmonised coloc data
# ==============================================================================

if (!file.exists(merged_path)) {
  stop("Missing ", merged_path, ". Run the peer-review coloc pipeline first.")
}

merged <- read_csv(merged_path, show_col_types = FALSE) %>%
  mutate(
    snp = as.character(snp),
    beta_exposure = as.numeric(beta_exposure),
    varbeta_exposure = as.numeric(varbeta_exposure),
    beta_outcome = as.numeric(beta_outcome),
    varbeta_outcome = as.numeric(varbeta_outcome)
  ) %>%
  filter(
    !is.na(snp),
    !is.na(beta_exposure),
    !is.na(varbeta_exposure),
    varbeta_exposure > 0,
    !is.na(beta_outcome),
    !is.na(varbeta_outcome),
    varbeta_outcome > 0
  ) %>%
  distinct(snp, .keep_all = TRUE)

message("Harmonised coloc SNPs available before LD lookup: ", nrow(merged))

# ==============================================================================
# 4. Get LD matrix from OpenGWAS reference panel
# ==============================================================================

message("Requesting EUR LD matrix from OpenGWAS reference panel...")

ld_raw <- tryCatch(
  {
    ieugwasr::ld_matrix(
      variants = merged$snp,
      pop = "EUR",
      with_alleles = FALSE
    )
  },
  error = function(e) {
    stop("LD matrix request failed: ", conditionMessage(e))
  }
)

ld_mat <- as.matrix(ld_raw)

# Defensive cleanup
storage.mode(ld_mat) <- "numeric"

common <- intersect(merged$snp, rownames(ld_mat))
common <- intersect(common, colnames(ld_mat))

message("SNPs retained in LD reference: ", length(common), " / ", nrow(merged))

if (length(common) < 10) {
  stop("Too few SNPs remained after LD lookup for SuSiE: ", length(common))
}

ld_mat <- ld_mat[common, common, drop = FALSE]

# Ensure exact order of data matches LD matrix
merged_ld <- merged %>%
  filter(snp %in% common) %>%
  mutate(order_index = match(snp, common)) %>%
  arrange(order_index) %>%
  select(-order_index)

ld_mat <- ld_mat[merged_ld$snp, merged_ld$snp, drop = FALSE]

# Clean matrix if needed
ld_mat[is.na(ld_mat)] <- 0
diag(ld_mat) <- 1

if (!all(rownames(ld_mat) == merged_ld$snp)) stop("LD row order does not match SNP order.")
if (!all(colnames(ld_mat) == merged_ld$snp)) stop("LD column order does not match SNP order.")

write_csv(as_tibble(ld_mat, rownames = "snp"), file.path(susie_dir, "gpnmb_coloc_susie_ld_matrix.csv"))
write_csv(merged_ld, file.path(results_dir, "gpnmb_coloc_susie_input_variants.csv"))

# ==============================================================================
# 5. Build SuSiE coloc datasets
# ==============================================================================

maf_complete <- "maf" %in% names(merged_ld) && all(!is.na(merged_ld$maf))

if (maf_complete) {
  message("Using MAF in SuSiE datasets.")
  D1 <- list(
    beta = merged_ld$beta_exposure,
    varbeta = merged_ld$varbeta_exposure,
    snp = merged_ld$snp,
    type = "quant",
    N = pqtl_N,
    MAF = merged_ld$maf,
    LD = ld_mat
  )
} else {
  message("MAF incomplete. Using sdY = 1 for pQTL SuSiE dataset.")
  D1 <- list(
    beta = merged_ld$beta_exposure,
    varbeta = merged_ld$varbeta_exposure,
    snp = merged_ld$snp,
    type = "quant",
    N = pqtl_N,
    sdY = 1,
    LD = ld_mat
  )
}

D2 <- list(
  beta = merged_ld$beta_outcome,
  varbeta = merged_ld$varbeta_outcome,
  snp = merged_ld$snp,
  type = "cc",
  N = pd_N,
  s = pd_s,
  LD = ld_mat
)

message("Checking coloc datasets with LD requirement...")
print(coloc::check_dataset(D1, req = "LD"))
print(coloc::check_dataset(D2, req = "LD"))

# ==============================================================================
# 6. Run SuSiE fine-mapping and coloc.susie
# ==============================================================================

message("Running SuSiE for pQTL exposure...")
S1 <- coloc::runsusie(D1)

message("Running SuSiE for PD outcome...")
S2 <- coloc::runsusie(D2)

exposure_summary <- summary(S1)
outcome_summary <- summary(S2)

print(exposure_summary)
print(outcome_summary)

# summary() may return a list-like object depending on package version.
# Capture safely.
try(write_csv(as_tibble(exposure_summary$cs), file.path(susie_dir, "gpnmb_coloc_susie_exposure_summary.csv")), silent = TRUE)
try(write_csv(as_tibble(outcome_summary$cs), file.path(susie_dir, "gpnmb_coloc_susie_outcome_summary.csv")), silent = TRUE)

message("Running coloc.susie...")
susie_coloc <- coloc::coloc.susie(S1, S2)

if (!is.null(susie_coloc$summary)) {
  susie_summary <- as_tibble(susie_coloc$summary)
} else {
  susie_summary <- as_tibble(susie_coloc)
}

print(susie_summary)
write_csv(susie_summary, file.path(results_dir, "gpnmb_coloc_susie_summary.csv"))

if (!is.null(susie_coloc$results)) {
  write_csv(as_tibble(susie_coloc$results), file.path(results_dir, "gpnmb_coloc_susie_results.csv"))
}

message("Done. Inspect gpnmb_coloc_susie_summary.csv.")
