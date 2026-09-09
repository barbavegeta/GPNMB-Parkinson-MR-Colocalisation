# GPNMB_coloc_peer_review_pipeline.R
#
# Peer-review oriented GPNMB coloc pipeline.
#
# What this script does:
#   1. Downloads dense regional pQTL and PD GWAS summary statistics from OpenGWAS.
#   2. Preserves chr, pos, alleles, beta, se, p, eaf/maf when available.
#   3. Deduplicates SNPs before joining.
#   4. Harmonises alleles and flips outcome beta when required.
#   5. Runs coloc.abf().
#   6. Runs a small prior-sensitivity analysis over p12 values.
#   7. Saves QC files, harmonised input, coloc summary, and SNP-level posterior results.
#
# Caveats:
#   - coloc.abf assumes one causal variant per trait in the analysed region.
#   - If MAF is missing, this script uses sdY = 1 for the quantitative pQTL trait.
#     That is only defensible if the protein phenotype is standardised.
#   - For a stronger peer-review standard, verify exposure scaling and run LD/multi-signal
#     colocalisation, e.g. coloc.susie, if LD is available.
#
# Required packages:
#   tidyverse, ieugwasr, coloc
#
# Run:
#   source("GPNMB_coloc_peer_review_pipeline.R")

library(tidyverse)
library(ieugwasr)
library(coloc)

# Allow running either from the package root or from 01_scripts.
if (basename(normalizePath(getwd(), winslash = "/", mustWork = FALSE)) == "01_scripts") {
  setwd("..")
}

results_dir <- file.path("02_results", "current_results")
regional_dir <- file.path(results_dir, "regional_inputs")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(regional_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 0. Parameters
# ==============================================================================

chr_locus <- "7"
start_pos <- 22814547
end_pos   <- 23814547

pqtl_id <- "prot-c-5080_131_3"
pd_id   <- "ieu-b-7"

pqtl_N <- 996

# Verify these against ieu-b-7 metadata before dissertation submission.
pd_ncase <- 33674
pd_ncontrol <- 449056
pd_N <- pd_ncase + pd_ncontrol
pd_s <- pd_ncase / pd_N

region_query <- paste0(chr_locus, ":", start_pos, "-", end_pos)

# ==============================================================================
# 1. Helper functions
# ==============================================================================

normalise_allele <- function(x) {
  toupper(as.character(x))
}

comp_allele <- function(x) {
  recode(
    normalise_allele(x),
    "A" = "T",
    "T" = "A",
    "C" = "G",
    "G" = "C",
    .default = NA_character_
  )
}

pick_col <- function(dat, candidates) {
  hit <- candidates[candidates %in% names(dat)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

optional_numeric <- function(dat, candidates) {
  hit <- pick_col(dat, candidates)
  if (is.na(hit)) return(rep(NA_real_, nrow(dat)))
  suppressWarnings(as.numeric(dat[[hit]]))
}

optional_character <- function(dat, candidates) {
  hit <- pick_col(dat, candidates)
  if (is.na(hit)) return(rep(NA_character_, nrow(dat)))
  as.character(dat[[hit]])
}

download_opengwas_region <- function(gwas_id, region, label) {
  message("Downloading ", label, " regional data from OpenGWAS: ", gwas_id)

  out <- tryCatch(
    {
      ieugwasr::associations(
        variants = region,
        id = gwas_id,
        proxies = 0
      )
    },
    error = function(e) {
      stop(label, " OpenGWAS download failed: ", conditionMessage(e))
    }
  )

  if (is.null(out) || nrow(out) == 0) {
    stop(label, " OpenGWAS query returned 0 rows.")
  }

  out <- as_tibble(out)
  message(label, " raw rows: ", nrow(out))
  message(label, " raw columns:")
  print(names(out))

  out
}

standardise_opengwas_region <- function(raw, label, N = NA_real_, ncase = NA_real_, ncontrol = NA_real_) {
  snp_col <- pick_col(raw, c("rsid", "snp", "SNP", "variant"))
  chr_col <- pick_col(raw, c("chr", "chromosome"))
  pos_col <- pick_col(raw, c("position", "pos", "bp"))
  beta_col <- pick_col(raw, c("beta", "b"))
  se_col <- pick_col(raw, c("se", "stderr", "standard_error"))
  p_col <- pick_col(raw, c("p", "pval", "pvalue"))
  ea_col <- pick_col(raw, c("ea", "effect_allele", "alt", "A1"))
  nea_col <- pick_col(raw, c("nea", "other_allele", "ref", "A2"))

  required <- c(snp_col, beta_col, se_col, p_col, ea_col, nea_col)
  if (any(is.na(required))) {
    detected <- tibble(
      field = c("snp", "beta", "se", "p", "effect_allele", "other_allele"),
      column = c(snp_col, beta_col, se_col, p_col, ea_col, nea_col)
    )
    print(detected)
    stop(label, " is missing one or more required OpenGWAS columns.")
  }

  eaf_vec <- optional_numeric(raw, c("eaf", "EAF", "af", "AF", "freq"))
  maf_vec <- optional_numeric(raw, c("maf", "MAF"))

  if (all(is.na(maf_vec)) && any(!is.na(eaf_vec))) {
    maf_vec <- pmin(eaf_vec, 1 - eaf_vec)
  }

  dat <- tibble(
    snp = as.character(raw[[snp_col]]),
    chr = optional_character(raw, c("chr", "chromosome")),
    pos = optional_numeric(raw, c("position", "pos", "bp")),
    effect_allele = normalise_allele(raw[[ea_col]]),
    other_allele = normalise_allele(raw[[nea_col]]),
    beta = as.numeric(raw[[beta_col]]),
    se = as.numeric(raw[[se_col]]),
    varbeta = as.numeric(raw[[se_col]])^2,
    p = as.numeric(raw[[p_col]]),
    eaf = eaf_vec,
    maf = maf_vec,
    n = N,
    ncase = ncase,
    ncontrol = ncontrol
  ) %>%
    filter(
      !is.na(snp),
      !is.na(beta),
      !is.na(se),
      !is.na(varbeta),
      varbeta > 0,
      !is.na(p),
      !is.na(effect_allele),
      !is.na(other_allele)
    )

  dup_n <- sum(duplicated(dat$snp))
  message(label, " duplicate SNP rows before deduplication: ", dup_n)

  dat <- dat %>%
    arrange(p, se) %>%
    distinct(snp, .keep_all = TRUE)

  message(label, " rows after cleaning/deduplication: ", nrow(dat))
  message(label, " rows with non-missing MAF: ", sum(!is.na(dat$maf)))

  dat
}

harmonise_regions <- function(exposure, outcome) {
  merged <- inner_join(
    exposure,
    outcome,
    by = "snp",
    suffix = c("_exposure", "_outcome")
  )

  message("Rows after join: ", nrow(merged))
  message("Unique SNPs after join: ", n_distinct(merged$snp))

  merged <- merged %>%
    arrange(p_exposure, p_outcome) %>%
    distinct(snp, .keep_all = TRUE)

  merged <- merged %>%
    mutate(
      same_alleles =
        effect_allele_exposure == effect_allele_outcome &
        other_allele_exposure == other_allele_outcome,

      reversed_alleles =
        effect_allele_exposure == other_allele_outcome &
        other_allele_exposure == effect_allele_outcome,

      same_complement =
        effect_allele_exposure == comp_allele(effect_allele_outcome) &
        other_allele_exposure == comp_allele(other_allele_outcome),

      reversed_complement =
        effect_allele_exposure == comp_allele(other_allele_outcome) &
        other_allele_exposure == comp_allele(effect_allele_outcome),

      allele_status = case_when(
        same_alleles ~ "same",
        reversed_alleles ~ "reversed",
        same_complement ~ "same_complement",
        reversed_complement ~ "reversed_complement",
        TRUE ~ "incompatible"
      ),

      beta_outcome_harmonised = case_when(
        allele_status %in% c("same", "same_complement") ~ beta_outcome,
        allele_status %in% c("reversed", "reversed_complement") ~ -beta_outcome,
        TRUE ~ NA_real_
      ),

      maf_for_coloc = case_when(
        !is.na(maf_exposure) ~ maf_exposure,
        !is.na(maf_outcome) ~ maf_outcome,
        TRUE ~ NA_real_
      )
    )

  message("Allele harmonisation status:")
  print(table(merged$allele_status, useNA = "ifany"))

  harmonised <- merged %>%
    filter(
      allele_status != "incompatible",
      !is.na(beta_outcome_harmonised)
    ) %>%
    transmute(
      snp = snp,
      chr = coalesce(chr_exposure, chr_outcome),
      pos = coalesce(pos_exposure, pos_outcome),

      beta_exposure = beta_exposure,
      se_exposure = se_exposure,
      varbeta_exposure = varbeta_exposure,
      p_exposure = p_exposure,
      effect_allele = effect_allele_exposure,
      other_allele = other_allele_exposure,

      beta_outcome = beta_outcome_harmonised,
      se_outcome = se_outcome,
      varbeta_outcome = varbeta_outcome,
      p_outcome = p_outcome,

      maf = maf_for_coloc,
      allele_status = allele_status
    ) %>%
    filter(
      !is.na(beta_exposure),
      !is.na(varbeta_exposure),
      !is.na(beta_outcome),
      !is.na(varbeta_outcome)
    ) %>%
    distinct(snp, .keep_all = TRUE)

  message("Rows retained after allele harmonisation: ", nrow(harmonised))
  message("Rows retained with MAF: ", sum(!is.na(harmonised$maf)))

  if (any(duplicated(harmonised$snp))) {
    stop("Duplicated SNPs remain after harmonisation.")
  }

  harmonised
}

run_coloc_with_prior <- function(dat, p12_value) {
  maf_complete <- all(!is.na(dat$maf))

  if (maf_complete) {
    d1 <- list(
      beta = dat$beta_exposure,
      varbeta = dat$varbeta_exposure,
      snp = dat$snp,
      type = "quant",
      N = pqtl_N,
      MAF = dat$maf
    )

    d2 <- list(
      beta = dat$beta_outcome,
      varbeta = dat$varbeta_outcome,
      snp = dat$snp,
      type = "cc",
      N = pd_N,
      s = pd_s,
      MAF = dat$maf
    )
  } else {
    warning("MAF incomplete. Using sdY = 1 for pQTL exposure for p12 = ", p12_value)

    d1 <- list(
      beta = dat$beta_exposure,
      varbeta = dat$varbeta_exposure,
      snp = dat$snp,
      type = "quant",
      N = pqtl_N,
      sdY = 1
    )

    d2 <- list(
      beta = dat$beta_outcome,
      varbeta = dat$varbeta_outcome,
      snp = dat$snp,
      type = "cc",
      N = pd_N,
      s = pd_s
    )
  }

  coloc::coloc.abf(
    dataset1 = d1,
    dataset2 = d2,
    p12 = p12_value
  )
}

# ==============================================================================
# 2. Download and standardise regional files
# ==============================================================================

pqtl_raw <- download_opengwas_region(pqtl_id, region_query, "pQTL")
pd_raw   <- download_opengwas_region(pd_id, region_query, "PD GWAS")

pqtl_region <- standardise_opengwas_region(
  pqtl_raw,
  label = "pQTL",
  N = pqtl_N
)

pd_region <- standardise_opengwas_region(
  pd_raw,
  label = "PD GWAS",
  N = pd_N,
  ncase = pd_ncase,
  ncontrol = pd_ncontrol
) %>%
  mutate(s = pd_s)

write_csv(pqtl_region, file.path(regional_dir, "gpnmb_pqtl_region_summary_peer_review.csv"))
write_csv(pd_region, file.path(regional_dir, "pd_gwas_gpnmb_region_summary_peer_review.csv"))

# ==============================================================================
# 3. Harmonise
# ==============================================================================

merged_harmonised <- harmonise_regions(pqtl_region, pd_region)

if (nrow(merged_harmonised) < 50) {
  warning("Fewer than 50 overlapping harmonised variants. Coloc may be underpowered or unstable.")
}

write_csv(merged_harmonised, file.path(results_dir, "gpnmb_coloc_merged_harmonised_peer_review.csv"))

# ==============================================================================
# 4. Main coloc and prior sensitivity
# ==============================================================================

p12_values <- c(1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4)

coloc_objects <- list()
sensitivity <- map_dfr(p12_values, function(p12v) {
  cr <- run_coloc_with_prior(merged_harmonised, p12v)
  coloc_objects[[as.character(p12v)]] <<- cr

  tibble(
    p12 = p12v,
    nsnps = unname(cr$summary[["nsnps"]]),
    PP.H0.abf = unname(cr$summary[["PP.H0.abf"]]),
    PP.H1.abf = unname(cr$summary[["PP.H1.abf"]]),
    PP.H2.abf = unname(cr$summary[["PP.H2.abf"]]),
    PP.H3.abf = unname(cr$summary[["PP.H3.abf"]]),
    PP.H4.abf = unname(cr$summary[["PP.H4.abf"]])
  )
})

write_csv(sensitivity, file.path(results_dir, "gpnmb_coloc_prior_sensitivity.csv"))
print(sensitivity)

# Use default coloc-like p12 = 1e-5 as the main result.
main_key <- "1e-05"
main_coloc <- coloc_objects[[main_key]]

coloc_summary_long <- tibble(
  Metric = names(main_coloc$summary),
  Value = as.numeric(main_coloc$summary)
)

coloc_summary_wide <- as_tibble(as.list(main_coloc$summary))

write_csv(coloc_summary_wide, file.path(results_dir, "gpnmb_final_coloc_results_peer_review.csv"))

message("Main coloc result:")
print(coloc_summary_wide)

message("Done. Check PP.H4.abf and prior sensitivity before interpreting.")
