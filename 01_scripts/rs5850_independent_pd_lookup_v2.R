library(TwoSampleMR)
library(ieugwasr)
library(readr)
library(dplyr)

# Allow running either from the package root or from 01_scripts.
if (basename(normalizePath(getwd(), winslash = "/", mustWork = FALSE)) == "01_scripts") {
  setwd("..")
}

# OpenGWAS token should be stored in ~/.Renviron as:
# OPENGWAS_JWT=your_token_here
cat("JWT present: ", nchar(ieugwasr::get_opengwas_jwt()) > 20, "\n")
print(ieugwasr::user())

outcomes <- c(
  "finn-b-G6_PARKINSON",
  "finn-b-PDSTRICT",
  "ebi-a-GCST90018894",
  "ebi-a-GCST90018674",
  "ieu-a-812",
  "ieu-a-818"
)

# Save candidate outcome metadata.
meta <- ieugwasr::gwasinfo(outcomes)
dir.create("03_supporting_inputs/independent_pd_outcomes", recursive = TRUE, showWarnings = FALSE)
write_csv(meta, "03_supporting_inputs/independent_pd_outcomes/opengwas_candidate_pd_outcome_metadata.csv")


# Proxy lookup where rs5850 is absent.
proxy <- extract_outcome_data(
  snps = "rs5850",
  outcomes = outcomes,
  proxies = TRUE,
  rsq = 0.8,
  align_alleles = 1,
  palindromes = 1,
  maf_threshold = 0.3
)

dir.create("02_results/current_results/somascan", recursive = TRUE, showWarnings = FALSE)
write_csv(proxy, "02_results/current_results/somascan/rs5850_independent_PD_proxy_lookup_v2.csv")

