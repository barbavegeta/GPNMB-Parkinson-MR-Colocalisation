options(stringsAsFactors = FALSE)
options(timeout = 600)

userlib <- path.expand("~/Rlibs")
dir.create(userlib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(userlib, .libPaths()))

pkgs <- c("readr", "dplyr", "ieugwasr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, lib = userlib)
  }
}

library(readr)
library(dplyr)
library(ieugwasr)

# ------------------------------------------------------------------
# Input/output paths for the final tidied dissertation package
# ------------------------------------------------------------------

# ============================================================
# Paths
# ============================================================

# RStudio may run this script while the working directory is 01_scripts.
# If so, move to the TO UPLOAD project root.
if (basename(normalizePath(getwd(), winslash = "/", mustWork = FALSE)) == "01_scripts") {
  setwd("..")
}

cat("Working directory:", getwd(), "\n\n")

input_file <- file.path(
  "02_results",
  "latest_olink",
  "gpnmb_ukbppp_ieu_b_7_harmonised_region_overlap.csv"
)

outdir <- file.path(
  "02_results",
  "latest_olink",
  "gpnmb_ld_clumped_mr_outputs"
)

if (!file.exists(input_file)) {
  stop(
    "Missing input file: ", input_file,
    "\nWorking directory: ", getwd()
  )
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("Input file:", input_file, "\n")
cat("Output directory:", outdir, "\n\n")

dat <- read_csv(input_file, show_col_types = FALSE)

required_cols <- c(
  "rsid", "ID", "ID_POS", "ALLELE0", "ALLELE1",
  "A1FREQ", "INFO", "BETA", "SE", "LOG10P",
  "out_beta_to_expA1", "out_se", "out_SS", "out_LP"
)

missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

dat2 <- dat %>%
  mutate(
    exposure_beta = BETA,
    exposure_se = SE,
    outcome_beta = out_beta_to_expA1,
    outcome_se = out_se,
    exposure_p = pmax(10^(-pmin(LOG10P, 300)), .Machine$double.xmin),
    F_stat = (exposure_beta / exposure_se)^2,
    outcome_N = out_SS
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
    exposure_p < 5e-8,
    outcome_N >= 450000
  ) %>%
  arrange(exposure_p)

write_csv(dat2, file.path(outdir, "01_filtered_genomewide_significant_harmonised_variants.csv"))

if (nrow(dat2) == 0) {
  stop("No variants survived QC + p < 5e-8 filter.")
}

run_clumped_mr <- function(dat_qc, clump_kb, clump_r2, label) {
  message("Running LD clumping: ", label)

  clump_input <- dat_qc %>%
    transmute(
      rsid = rsid,
      pval = exposure_p,
      id = "GPNMB_UKBPPP_Olink"
    ) %>%
    distinct(rsid, .keep_all = TRUE)

  clumped_ids <- ieugwasr::ld_clump(
    dat = clump_input,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    clump_p = 5e-8,
    pop = "EUR"
  )

  clumped <- dat_qc %>%
    semi_join(clumped_ids, by = "rsid") %>%
    arrange(exposure_p)

  write_csv(
    clumped,
    file.path(outdir, paste0("02_clumped_instruments_", label, ".csv"))
  )

  if (nrow(clumped) == 0) {
    res <- tibble(
      setting = label,
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      n_instruments = 0,
      mr_method = "No instruments retained after LD clumping",
      mr_beta = NA_real_,
      mr_se = NA_real_,
      mr_p = NA_real_,
      mr_OR = NA_real_,
      mr_OR_lci = NA_real_,
      mr_OR_uci = NA_real_
    )
    return(res)
  }

  if (nrow(clumped) == 1) {
    d <- clumped[1, ]

    beta <- d$outcome_beta / d$exposure_beta
    se <- sqrt(
      (d$outcome_se^2 / d$exposure_beta^2) +
        ((d$outcome_beta^2 * d$exposure_se^2) / d$exposure_beta^4)
    )

    res <- tibble(
      setting = label,
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      n_instruments = 1,
      instruments = d$rsid,
      lead_variant = d$rsid,
      mr_method = "LD-clumped single-instrument Wald ratio",
      mr_beta = beta,
      mr_se = se,
      mr_p = 2 * pnorm(abs(beta / se), lower.tail = FALSE),
      mr_OR = exp(beta),
      mr_OR_lci = exp(beta - 1.96 * se),
      mr_OR_uci = exp(beta + 1.96 * se)
    )

    return(res)
  }

  bx <- clumped$exposure_beta
  by <- clumped$outcome_beta
  sey <- clumped$outcome_se
  w <- 1 / sey^2

  beta_ivw <- sum(w * bx * by) / sum(w * bx^2)
  se_ivw <- sqrt(1 / sum(w * bx^2))
  z_ivw <- beta_ivw / se_ivw
  p_ivw <- 2 * pnorm(abs(z_ivw), lower.tail = FALSE)

  q <- sum(w * (by - beta_ivw * bx)^2)
  q_df <- nrow(clumped) - 1
  q_p <- pchisq(q, df = q_df, lower.tail = FALSE)
  
  # Multiplicative random-effects IVW
  #
  # The IVW point estimate remains unchanged, but its standard error is
  # multiplied by the residual standard deviation sqrt(Q / df).
  # max(1, ...) prevents under-dispersion from making the random-effects
  # estimate artificially more precise than fixed-effect IVW.
  
  re_phi <- max(1, q / q_df)
  re_se <- se_ivw * sqrt(re_phi)
  re_z <- beta_ivw / re_se
  re_p <- 2 * pnorm(abs(re_z), lower.tail = FALSE)
  
  re_OR <- exp(beta_ivw)
  re_OR_lci <- exp(beta_ivw - 1.96 * re_se)
  re_OR_uci <- exp(beta_ivw + 1.96 * re_se)

  ratio <- by / bx
  ratio_se <- sqrt(
    (clumped$outcome_se^2 / clumped$exposure_beta^2) +
      ((clumped$outcome_beta^2 * clumped$exposure_se^2) / clumped$exposure_beta^4)
  )

  per_variant <- clumped %>%
    mutate(
      ratio_beta = ratio,
      ratio_se = ratio_se,
      ratio_p = 2 * pnorm(abs(ratio_beta / ratio_se), lower.tail = FALSE),
      ratio_OR = exp(ratio_beta),
      ratio_OR_lci = exp(ratio_beta - 1.96 * ratio_se),
      ratio_OR_uci = exp(ratio_beta + 1.96 * ratio_se)
    )

  if (label != "sensitivity_r2_0p01_1Mb") {
    write_csv(
      per_variant,
      file.path(outdir, paste0("03_per_variant_ratios_", label, ".csv"))
    )
  }

  res <- tibble(
    setting = label,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    n_instruments = nrow(clumped),
    instruments = paste(clumped$rsid, collapse = ";"),
    lead_variant = clumped$rsid[1],
    mr_method = "LD-clumped fixed-effect IVW",
    mr_beta = beta_ivw,
    mr_se = se_ivw,
    mr_p = p_ivw,
    mr_OR = exp(beta_ivw),
    mr_OR_lci = exp(beta_ivw - 1.96 * se_ivw),
    mr_OR_uci = exp(beta_ivw + 1.96 * se_ivw),
    cochran_Q = q,
    cochran_Q_df = q_df,
    cochran_Q_p = q_p,
    
    multiplicative_RE_phi = re_phi,
    multiplicative_RE_beta = beta_ivw,
    multiplicative_RE_se = re_se,
    multiplicative_RE_p = re_p,
    multiplicative_RE_OR = re_OR,
    multiplicative_RE_OR_lci = re_OR_lci,
    multiplicative_RE_OR_uci = re_OR_uci,
    
    median_ratio_beta = median(ratio, na.rm = TRUE),
    median_ratio_OR = exp(median(ratio, na.rm = TRUE))
  )

  return(res)
}

results <- bind_rows(
  run_clumped_mr(dat2, clump_kb = 10000, clump_r2 = 0.001, label = "strict_r2_0p001_10Mb"),
  run_clumped_mr(dat2, clump_kb = 10000, clump_r2 = 0.01,  label = "standard_r2_0p01_10Mb"),
  run_clumped_mr(dat2, clump_kb = 1000,  clump_r2 = 0.01,  label = "sensitivity_r2_0p01_1Mb"),
  run_clumped_mr(dat2, clump_kb = 1000,  clump_r2 = 0.05,  label = "liberal_r2_0p05_1Mb")
)

write_csv(results, file.path(outdir, "04_all_ld_clumped_mr_results.csv"))

random_effects_results <- results %>%
  filter(n_instruments >= 2) %>%
  select(
    setting,
    clump_kb,
    clump_r2,
    n_instruments,
    instruments,
    lead_variant,
    cochran_Q,
    cochran_Q_df,
    cochran_Q_p,
    multiplicative_RE_phi,
    multiplicative_RE_beta,
    multiplicative_RE_se,
    multiplicative_RE_p,
    multiplicative_RE_OR,
    multiplicative_RE_OR_lci,
    multiplicative_RE_OR_uci
  )

write_csv(
  random_effects_results,
  file.path(outdir, "05_multiplicative_random_effects_ivw_results.csv")
)

sink(file.path(outdir, "00_READ_ME_RESULTS_SUMMARY.txt"))
cat("GPNMB corrected Olink / UKB-PPP -> ieu-b-7 LD-clumped MR\n")
cat("Created:", as.character(Sys.time()), "\n\n")
cat("Input file:", input_file, "\n")
cat("Raw harmonised variants:", nrow(dat), "\n")
cat("After QC + p < 5e-8:", nrow(dat2), "\n\n")
cat("MR results:\n")
print(results)
cat("\nInterpretation rules:\n")
cat("- If n_instruments = 1, that setting is a clumped single-instrument Wald ratio, not IVW.\n")
cat("- If n_instruments >= 2, that setting is LD-clumped fixed-effect IVW.\n")
cat("- Use strict or standard clumping as the main analysis; liberal clumping is sensitivity only.\n")
cat("- Do not use all regional variants as independent instruments.\n")
cat("\nMultiplicative random-effects IVW results:\n")

print(
  results %>%
    filter(n_instruments >= 2) %>%
    select(
      setting,
      n_instruments,
      multiplicative_RE_beta,
      multiplicative_RE_se,
      multiplicative_RE_p,
      multiplicative_RE_OR,
      multiplicative_RE_OR_lci,
      multiplicative_RE_OR_uci,
      cochran_Q,
      cochran_Q_df,
      cochran_Q_p
    )
)

cat("\nRandom-effects interpretation:\n")
cat("- Multiplicative random-effects IVW retains the fixed-effect IVW point estimate.\n")
cat("- The standard error is inflated by sqrt(max(1, Cochran Q / Q df)).\n")
cat("- This accounts for over-dispersion between SNP-specific estimates but does not explain the source of heterogeneity.\n")
sink()

cat("\nDONE. Outputs written to:\n")
cat(outdir, "\n")
print(list.files(outdir, full.names = TRUE))

