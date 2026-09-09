# GPNMB corrected Olink / UKB-PPP -> ieu-b-7 Parkinson's MR
# Input needed: gpnmb_processed_outputs.zip or extracted gpnmb_ukbppp_ieu_b_7_harmonised_region_overlap.csv

options(stringsAsFactors = FALSE)

needed <- c("readr", "dplyr")
for (p in needed) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
library(readr)
library(dplyr)

# Allow running either from the package root or from 01_scripts.
if (basename(normalizePath(getwd(), winslash = "/", mustWork = FALSE)) == "01_scripts") {
  setwd("..")
}

outdir <- file.path("02_results", "latest_olink", "gpnmb_corrected_mr_outputs")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 1. Load final harmonised regional input
input_file <- file.path("02_results", "latest_olink", "gpnmb_ukbppp_ieu_b_7_harmonised_region_overlap.csv")
if (!file.exists(input_file)) stop("Missing input file: ", input_file)
dat <- read_csv(input_file, show_col_types = FALSE)

# 2. Clean and filter
dat2 <- dat %>%
  mutate(
    exposure_beta = BETA,
    exposure_se   = SE,
    outcome_beta  = out_beta_to_expA1,
    outcome_se    = out_se,
    exposure_p    = 10^(-pmin(LOG10P, 300)),
    F_stat        = (exposure_beta / exposure_se)^2,
    outcome_N     = out_SS,
    outcome_cases = out_NC
  ) %>%
  filter(
    !is.na(rsid),
    !is.na(exposure_beta),
    !is.na(exposure_se),
    !is.na(outcome_beta),
    !is.na(outcome_se),
    exposure_se > 0,
    outcome_se > 0,
    INFO >= 0.8,
    A1FREQ >= 0.01,
    A1FREQ <= 0.99,
    F_stat >= 10,
    outcome_N >= 450000
  ) %>%
  arrange(desc(LOG10P))

# 3. Lead-variant MR: strongest exposure variant after QC
lead <- dat2 %>% slice(1)

wald_mr <- function(d) {
  beta <- d$outcome_beta / d$exposure_beta

  # Delta-method SE, including exposure uncertainty
  se <- sqrt(
    (d$outcome_se^2 / d$exposure_beta^2) +
      ((d$outcome_beta^2 * d$exposure_se^2) / d$exposure_beta^4)
  )

  tibble(
    rsid = d$rsid,
    ID = d$ID,
    ID_POS = d$ID_POS,
    effect_allele = d$ALLELE1,
    other_allele = d$ALLELE0,
    A1FREQ = d$A1FREQ,
    INFO = d$INFO,
    exposure_beta = d$exposure_beta,
    exposure_se = d$exposure_se,
    exposure_LOG10P = d$LOG10P,
    F_stat = d$F_stat,
    outcome_beta = d$outcome_beta,
    outcome_se = d$outcome_se,
    outcome_LP = d$out_LP,
    outcome_N = d$outcome_N,
    mr_method = "Lead-variant Wald ratio",
    mr_beta = beta,
    mr_se = se,
    mr_p = 2 * pnorm(abs(beta / se), lower.tail = FALSE),
    mr_OR = exp(beta),
    mr_OR_lci = exp(beta - 1.96 * se),
    mr_OR_uci = exp(beta + 1.96 * se)
  )
}

lead_res <- wald_mr(lead)
write_csv(lead_res, file.path(outdir, "01_lead_variant_wald_mr.csv"))

cat("Lead-variant MR complete.\n")
