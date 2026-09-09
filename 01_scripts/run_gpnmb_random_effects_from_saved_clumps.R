options(stringsAsFactors = FALSE)

userlib <- path.expand("~/Rlibs")
dir.create(userlib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(userlib, .libPaths()))

pkgs <- c("readr", "dplyr")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, lib = userlib, repos = "https://cloud.r-project.org")
  }
}

library(readr)
library(dplyr)

# Allow running either from the package root or from 01_scripts.
if (basename(normalizePath(getwd(), winslash = "/", mustWork = FALSE)) == "01_scripts") {
  setwd("..")
}

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

outdir <- file.path(
  "02_results",
  "latest_olink",
  "gpnmb_ld_clumped_mr_outputs"
)

if (!dir.exists(outdir)) {
  stop("Missing output directory: ", outdir)
}

# ------------------------------------------------------------
# Calculate FE-IVW and multiplicative random-effects IVW
# directly from EXISTING final LD-clumped instruments.
#
# No OpenGWAS call is performed here.
# ------------------------------------------------------------

calculate_mr <- function(filename, setting, clump_kb, clump_r2) {

  f <- file.path(outdir, filename)

  if (!file.exists(f)) {
    stop("Missing clumped instrument file: ", f)
  }

  dat <- read_csv(f, show_col_types = FALSE)

  required <- c(
    "rsid",
    "exposure_beta",
    "outcome_beta",
    "outcome_se"
  )

  missing <- setdiff(required, names(dat))

  if (length(missing) > 0) {
    stop(
      "Missing columns in ", filename, ": ",
      paste(missing, collapse = ", ")
    )
  }

  if (nrow(dat) < 2) {
    stop(
      "Cannot calculate IVW for ", setting,
      ": fewer than two instruments."
    )
  }

  bx <- dat$exposure_beta
  by <- dat$outcome_beta
  sey <- dat$outcome_se

  # Fixed-effect IVW
  w <- 1 / sey^2

  beta_ivw <- sum(w * bx * by) / sum(w * bx^2)
  se_ivw <- sqrt(1 / sum(w * bx^2))
  z_ivw <- beta_ivw / se_ivw
  p_ivw <- 2 * pnorm(abs(z_ivw), lower.tail = FALSE)

  # Cochran Q
  q <- sum(w * (by - beta_ivw * bx)^2)
  q_df <- nrow(dat) - 1
  q_p <- pchisq(q, df = q_df, lower.tail = FALSE)

  # Multiplicative random-effects IVW
  phi <- max(1, q / q_df)
  re_se <- se_ivw * sqrt(phi)
  re_z <- beta_ivw / re_se
  re_p <- 2 * pnorm(abs(re_z), lower.tail = FALSE)

  tibble(
    setting = setting,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    n_instruments = nrow(dat),
    instruments = paste(dat$rsid, collapse = ";"),

    fixed_effect_beta = beta_ivw,
    fixed_effect_se = se_ivw,
    fixed_effect_p = p_ivw,
    fixed_effect_OR = exp(beta_ivw),
    fixed_effect_OR_lci = exp(beta_ivw - 1.96 * se_ivw),
    fixed_effect_OR_uci = exp(beta_ivw + 1.96 * se_ivw),

    cochran_Q = q,
    cochran_Q_df = q_df,
    cochran_Q_p = q_p,

    multiplicative_RE_phi = phi,
    multiplicative_RE_beta = beta_ivw,
    multiplicative_RE_se = re_se,
    multiplicative_RE_p = re_p,
    multiplicative_RE_OR = exp(beta_ivw),
    multiplicative_RE_OR_lci = exp(beta_ivw - 1.96 * re_se),
    multiplicative_RE_OR_uci = exp(beta_ivw + 1.96 * re_se)
  )
}

# ------------------------------------------------------------
# Use the four SAVED final clumping solutions
# ------------------------------------------------------------

results <- bind_rows(

  calculate_mr(
    "02_clumped_instruments_strict_r2_0p001_10Mb.csv",
    "strict_r2_0p001_10Mb",
    10000,
    0.001
  ),

  calculate_mr(
    "02_clumped_instruments_standard_r2_0p01_10Mb.csv",
    "standard_r2_0p01_10Mb",
    10000,
    0.01
  ),

  calculate_mr(
    "02_clumped_instruments_sensitivity_r2_0p01_1Mb.csv",
    "sensitivity_r2_0p01_1Mb",
    1000,
    0.01
  ),

  calculate_mr(
    "02_clumped_instruments_liberal_r2_0p05_1Mb.csv",
    "liberal_r2_0p05_1Mb",
    1000,
    0.05
  )
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

outfile <- file.path(
  outdir,
  "05_multiplicative_random_effects_ivw_results.csv"
)

write_csv(results, outfile)

cat("\n============================================\n")
cat("Multiplicative random-effects IVW completed\n")
cat("============================================\n\n")

print(
  results %>%
    select(
      setting,
      n_instruments,
      fixed_effect_beta,
      fixed_effect_se,
      fixed_effect_p,
      cochran_Q,
      cochran_Q_df,
      cochran_Q_p,
      multiplicative_RE_se,
      multiplicative_RE_p,
      multiplicative_RE_OR,
      multiplicative_RE_OR_lci,
      multiplicative_RE_OR_uci
    ),
  n = Inf,
  width = Inf
)

cat("\nWritten to:\n", outfile, "\n")
