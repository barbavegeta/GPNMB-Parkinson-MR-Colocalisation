# GPNMB_MR_full_pipeline_with_olink_rsid_mapping.R
#
# Full pipeline:
#   1. Format Olink GPNMB credible-set cis-pQTL variants.
#   2. Reconstruct Olink SE from beta and p-value because the credible-set export has no SE column.
#   3. Load already-formatted SomaScan GPNMB exposure data safely, forcing allele columns to character.
#   4. Extract Parkinson's disease outcome data from OpenGWAS.
#   5. Run SomaScan single-SNP MR.
#   6. Try to map Olink chr_pos_ref_alt variant IDs to rsIDs via Ensembl Variant Recoder.
#   7. Retry Olink outcome extraction and MR using mapped rsIDs.
#
# Run this from the folder containing:
#   - ENSG00000136235-qtl-credible-sets-target.tsv
#   - gpnmb_exposure_prot-c-5080_131_3.csv
#
# Required packages:
#   tidyverse, TwoSampleMR, httr, jsonlite
#
# Important caveat:
#   The Olink credible-set export uses Open Targets GRCh38 chr_pos_ref_alt IDs.
#   OpenGWAS usually matches outcome variants more reliably with rsIDs.
#   The Ensembl mapping step is therefore a best-effort rescue, not guaranteed to map every variant.

library(tidyverse)
library(TwoSampleMR)
library(httr)
library(jsonlite)

# Allow running either from the package root or from 01_scripts.
if (basename(normalizePath(getwd(), winslash = "/", mustWork = FALSE)) == "01_scripts") {
  setwd("..")
}

results_dir <- file.path("02_results", "current_results")
olink_dir <- file.path(results_dir, "olink")
somascan_dir <- file.path(results_dir, "somascan")
dir.create(olink_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(somascan_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 0. Helper functions
# ==============================================================================

clean_pvalue <- function(x) {
  x_chr <- x %>%
    as.character() %>%
    str_replace_all("−", "-") %>%
    str_replace_all("<", "") %>%
    str_replace_all(">", "") %>%
    str_replace_all("\\s+", "") %>%
    str_replace_all("([0-9.]+)[xX×]10\\^?([+-]?[0-9]+)", "\\1e\\2")

  suppressWarnings(as.numeric(x_chr))
}

stop_if_empty <- function(dat, label) {
  if (is.null(dat) || nrow(dat) == 0) {
    stop(label, " has 0 rows. Stop here.")
  }
}

safe_nrow <- function(dat) {
  if (is.null(dat)) 0 else nrow(dat)
}

extract_outcome_safely <- function(snps, outcome_id, label, proxies = FALSE) {
  if (length(snps) == 0) {
    message(label, ": no SNPs supplied to extract_outcome_data().")
    return(NULL)
  }

  out <- tryCatch(
    {
      extract_outcome_data(
        snps = unique(snps),
        outcomes = outcome_id,
        proxies = proxies
      )
    },
    error = function(e) {
      message(label, ": outcome extraction failed: ", conditionMessage(e))
      return(NULL)
    }
  )

  out
}

run_mr_safely <- function(harmonised_dat, label, harmonised_path, odds_path) {
  if (is.null(harmonised_dat) || nrow(harmonised_dat) == 0) {
    message(label, ": harmonised dataset has 0 rows. MR skipped.")
    return(NULL)
  }

  if (!"mr_keep" %in% names(harmonised_dat)) {
    message(label, ": no mr_keep column. MR skipped.")
    return(NULL)
  }

  usable <- harmonised_dat %>% filter(mr_keep == TRUE)

  if (nrow(usable) == 0) {
    message(label, ": 0 SNPs retained after harmonisation. MR skipped.")
    return(NULL)
  }

  write_csv(usable, harmonised_path)

  mr_res <- tryCatch(
    {
      mr(usable)
    },
    error = function(e) {
      message(label, ": MR failed: ", conditionMessage(e))
      return(NULL)
    }
  )

  if (is.null(mr_res) || nrow(mr_res) == 0 || !"b" %in% names(mr_res)) {
    message(label, ": MR returned no valid estimates. Odds ratios skipped.")
    if (!is.null(mr_res) && is.data.frame(mr_res)) {
    }
    return(list(harmonised = usable, mr = mr_res, or = NULL))
  }

  or_res <- generate_odds_ratios(mr_res)

  write_csv(or_res, odds_path)

  cat("\n", label, " MR results:\n", sep = "")
  print(mr_res)
  cat("\n", label, " odds ratios:\n", sep = "")
  print(or_res)

  list(harmonised = usable, mr = mr_res, or = or_res)
}

summarise_mr_result <- function(result_object, platform, exposure_label, outcome_label, outcome_id) {
  if (is.null(result_object) ||
      is.null(result_object$mr) ||
      is.null(result_object$or) ||
      nrow(result_object$mr) == 0 ||
      nrow(result_object$or) == 0) {
    return(tibble())
  }

  tibble(
    platform = platform,
    exposure = exposure_label,
    outcome = outcome_label,
    outcome_id = outcome_id,
    method = result_object$mr$method,
    nsnp = result_object$mr$nsnp,
    beta = result_object$mr$b,
    se = result_object$mr$se,
    pval = result_object$mr$pval,
    OR = result_object$or$or,
    OR_lci95 = result_object$or$or_lci95,
    OR_uci95 = result_object$or$or_uci95
  )
}

# Ensembl Variant Recoder mapping helper.
# This tries several possible representations and pulls rsIDs from the JSON response.
query_ensembl_variant_recoder <- function(query_id) {
  url <- paste0(
    "https://rest.ensembl.org/variant_recoder/homo_sapiens/",
    URLencode(query_id, reserved = TRUE),
    "?content-type=application/json"
  )

  res <- tryCatch(
    {
      GET(
        url,
        add_headers(
          "Accept" = "application/json",
          "Content-Type" = "application/json"
        ),
        timeout(20)
      )
    },
    error = function(e) {
      return(NULL)
    }
  )

  if (is.null(res) || status_code(res) != 200) {
    return(character())
  }

  txt <- content(res, as = "text", encoding = "UTF-8")

  rsids <- str_extract_all(txt, "rs[0-9]+")[[1]] %>%
    unique()

  rsids
}

map_one_opentargets_variant_to_rsid <- function(chr, pos, ref, alt, original_id) {
  chr <- as.character(chr)
  pos <- as.integer(pos)
  ref <- as.character(ref)
  alt <- as.character(alt)

  queries <- c(
    original_id,
    str_replace_all(original_id, "_", "-")
  )

  # HGVS genomic syntax is simple for SNVs.
  # Indel HGVS notation is more complex, so for indels we still try original/dash/SPDI-style queries.
  if (nchar(ref) == 1 && nchar(alt) == 1) {
    queries <- c(
      queries,
      paste0(chr, ":g.", pos, ref, ">", alt)
    )
  }

  # GRCh38 RefSeq accession for chromosome 7.
  # SPDI is 0-based, but services vary in accepted forms; try both pos and pos-1.
  if (chr == "7") {
    queries <- c(
      queries,
      paste0("NC_000007.14:", pos, ":", ref, ":", alt),
      paste0("NC_000007.14:", pos - 1, ":", ref, ":", alt)
    )
  }

  queries <- unique(queries)

  for (q in queries) {
    Sys.sleep(0.2)
    rsids <- query_ensembl_variant_recoder(q)

    if (length(rsids) > 0) {
      return(tibble(
        SNP_original = original_id,
        chr = chr,
        pos = pos,
        ref = ref,
        alt = alt,
        query_used = q,
        rsid = rsids[1],
        all_rsids = paste(rsids, collapse = ";")
      ))
    }
  }

  tibble(
    SNP_original = original_id,
    chr = chr,
    pos = pos,
    ref = ref,
    alt = alt,
    query_used = NA_character_,
    rsid = NA_character_,
    all_rsids = NA_character_
  )
}

# ==============================================================================
# 1. Olink credible-set exposure formatting
# ==============================================================================

credible_sets <- read_tsv(
  file.path("03_supporting_inputs", "OpenTargets", "ENSG00000136235-qtl-credible-sets-target.tsv"),
  show_col_types = FALSE
)

cat("Credible-set rows:", nrow(credible_sets), "\n")
cat("Credible-set columns:\n")
print(names(credible_sets))

olink_raw <- credible_sets %>%
  filter(study == "UKB_PPP_EUR_GPNMB_Q14956_OID20173_v1") %>%
  filter(str_detect(type, "pQTL")) %>%
  filter(str_detect(type, "isTrans:false"))

cat("Olink cis-pQTL rows:", nrow(olink_raw), "\n")
stop_if_empty(olink_raw, "olink_raw")

olink_parsed <- olink_raw %>%
  separate(
    leadVariant,
    into = c("chr", "pos", "ref_allele", "alt_allele"),
    sep = "_",
    remove = FALSE,
    convert = FALSE
  ) %>%
  mutate(
    chr = as.character(chr),
    pos = as.integer(pos),
    beta = as.numeric(beta),
    effect_allele = as.character(alt_allele),
    other_allele = as.character(ref_allele),
    pValue_original = pValue,
    pValue = clean_pvalue(pValue),
    sample_size = 34151,
    exposure_name = "GPNMB (Olink UKB-PPP)"
  )

bad_p <- olink_parsed %>%
  filter(is.na(pValue) | pValue <= 0) %>%
  select(leadVariant, pValue_original, pValue)

if (nrow(bad_p) > 0) {
  print(bad_p, n = Inf)
  stop("Some Olink p-values could not be parsed.")
}

# Reconstruct SE because the Olink credible-set export lacks an SE column.
olink_parsed <- olink_parsed %>%
  mutate(
    pValue_safe = pmax(pValue, .Machine$double.xmin),
    z = qnorm(pValue_safe / 2, lower.tail = FALSE),
    se = abs(beta) / z
  ) %>%
  filter(
    !is.na(beta),
    !is.na(se),
    is.finite(se),
    se > 0,
    !is.na(effect_allele),
    !is.na(other_allele)
  )

cat("Olink usable rows after SE reconstruction:", nrow(olink_parsed), "\n")
stop_if_empty(olink_parsed, "olink_parsed")

olink_exposure_dat <- format_data(
  olink_parsed,
  type = "exposure",
  snp_col = "leadVariant",
  chr_col = "chr",
  pos_col = "pos",
  beta_col = "beta",
  se_col = "se",
  pval_col = "pValue",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  samplesize_col = "sample_size",
  phenotype_col = "exposure_name"
)

cat("Formatted Olink exposure rows:", nrow(olink_exposure_dat), "\n")

write_csv(olink_parsed, file.path(olink_dir, "olink_gpnmb_cis_pqtl_reconstructed_se.csv"))

# ==============================================================================
# 2. SomaScan exposure loading
# ==============================================================================

# This CSV is already in TwoSampleMR format.
# Force allele columns to character so T is not read as TRUE.
somascan_exposure_dat <- read_csv(
  file.path(somascan_dir, "gpnmb_exposure_prot-c-5080_131_3.csv"),
  show_col_types = FALSE,
  col_types = cols(
    id.exposure = col_character(),
    chr.exposure = col_character(),
    pos.exposure = col_double(),
    SNP = col_character(),
    effect_allele.exposure = col_character(),
    other_allele.exposure = col_character(),
    eaf.exposure = col_double(),
    beta.exposure = col_double(),
    se.exposure = col_double(),
    pval.exposure = col_double(),
    samplesize.exposure = col_double(),
    exposure = col_character(),
    mr_keep.exposure = col_logical(),
    pval_origin.exposure = col_character(),
    data_source.exposure = col_character()
  )
)

cat("SomaScan exposure rows:", nrow(somascan_exposure_dat), "\n")
stop_if_empty(somascan_exposure_dat, "somascan_exposure_dat")

cat("SomaScan allele check:\n")
somascan_exposure_dat %>%
  select(SNP, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure, pval.exposure) %>%
  print()


# ==============================================================================
# 3. Outcome extraction: original Olink IDs and SomaScan rsID
# ==============================================================================

pd_outcome_id <- "ieu-b-7"
pd_outcome_label <- "Parkinson's disease"

outcome_olink_chrpos <- extract_outcome_safely(
  snps = unique(olink_exposure_dat$SNP),
  outcome_id = pd_outcome_id,
  label = "Olink chr_pos_ref_alt",
  proxies = FALSE
)

outcome_somascan <- extract_outcome_safely(
  snps = unique(somascan_exposure_dat$SNP),
  outcome_id = pd_outcome_id,
  label = "SomaScan",
  proxies = FALSE
)

cat("Olink chr_pos_ref_alt outcome rows:", safe_nrow(outcome_olink_chrpos), "\n")
cat("SomaScan outcome rows:", safe_nrow(outcome_somascan), "\n")


# ==============================================================================
# 4. Harmonise and run SomaScan MR
# ==============================================================================

somascan_harmonised <- NULL
somascan_results <- NULL

if (!is.null(outcome_somascan) && nrow(outcome_somascan) > 0) {
  somascan_harmonised <- harmonise_data(
    exposure_dat = somascan_exposure_dat,
    outcome_dat = outcome_somascan,
    action = 2
  )

  cat("SomaScan harmonised rows:", nrow(somascan_harmonised), "\n")
  print(table(somascan_harmonised$mr_keep))

  cat("SomaScan harmonisation check:\n")
  somascan_harmonised %>%
    select(
      SNP,
      effect_allele.exposure,
      other_allele.exposure,
      effect_allele.outcome,
      other_allele.outcome,
      beta.exposure,
      beta.outcome,
      mr_keep,
      remove,
      palindromic,
      ambiguous
    ) %>%
    print()

  somascan_results <- run_mr_safely(
    harmonised_dat = somascan_harmonised,
    label = "SomaScan",
    harmonised_path = file.path(somascan_dir, "somascan_gpnmb_pd_harmonised.csv"),
    odds_path = file.path(somascan_dir, "somascan_gpnmb_pd_odds_ratios.csv")
  )
} else {
  message("SomaScan skipped: no Parkinson's outcome rows found.")
}

# ==============================================================================
# 5. Olink rsID mapping via Ensembl Variant Recoder
# ==============================================================================

olink_variant_input <- olink_exposure_dat %>%
  distinct(SNP) %>%
  separate(
    SNP,
    into = c("chr", "pos", "ref", "alt"),
    sep = "_",
    remove = FALSE,
    convert = FALSE
  ) %>%
  mutate(pos = as.integer(pos))

cat("Attempting Ensembl rsID mapping for Olink variants:", nrow(olink_variant_input), "variants\n")

olink_rsid_map <- pmap_dfr(
  list(
    chr = olink_variant_input$chr,
    pos = olink_variant_input$pos,
    ref = olink_variant_input$ref,
    alt = olink_variant_input$alt,
    original_id = olink_variant_input$SNP
  ),
  map_one_opentargets_variant_to_rsid
)

write_csv(olink_rsid_map, file.path(olink_dir, "olink_variant_to_rsid_map.csv"))

cat("Olink rsID mapping summary:\n")
olink_rsid_map %>%
  summarise(
    total_variants = n(),
    mapped_to_rsid = sum(!is.na(rsid)),
    unmapped = sum(is.na(rsid))
  ) %>%
  print()

cat("Olink rsID map:\n")
olink_rsid_map %>%
  as_tibble() %>%
  print(n = Inf)

# ==============================================================================
# 6. Retry Olink MR with mapped rsIDs
# ==============================================================================

olink_exposure_rsid_dat <- olink_exposure_dat %>%
  rename(SNP_original = SNP) %>%
  left_join(
    olink_rsid_map %>% select(SNP_original, rsid),
    by = "SNP_original"
  ) %>%
  filter(!is.na(rsid)) %>%
  mutate(SNP = rsid) %>%
  select(-rsid) %>%
  relocate(SNP, .before = SNP_original)

if (nrow(olink_exposure_rsid_dat) > 0) {
  # If several Open Targets variants map to the same rsID, keep the strongest p-value.
  olink_exposure_rsid_dat <- olink_exposure_rsid_dat %>%
    arrange(pval.exposure) %>%
    distinct(SNP, .keep_all = TRUE)

  write_csv(olink_exposure_rsid_dat, file.path(olink_dir, "olink_gpnmb_exposure_twosamplemr_rsid_mapped.csv"))
}

cat("Olink exposure rows after rsID mapping:", nrow(olink_exposure_rsid_dat), "\n")

outcome_olink_rsid <- NULL
olink_harmonised_rsid <- NULL
olink_results_rsid <- NULL

if (nrow(olink_exposure_rsid_dat) > 0) {
  outcome_olink_rsid <- extract_outcome_safely(
    snps = unique(olink_exposure_rsid_dat$SNP),
    outcome_id = pd_outcome_id,
    label = "Olink rsID-mapped",
    proxies = FALSE
  )

  cat("Olink rsID-mapped outcome rows:", safe_nrow(outcome_olink_rsid), "\n")

  if (!is.null(outcome_olink_rsid) && nrow(outcome_olink_rsid) > 0) {
    write_csv(outcome_olink_rsid, file.path(olink_dir, "olink_gpnmb_pd_outcome_rows_rsid_mapped.csv"))

    olink_harmonised_rsid <- harmonise_data(
      exposure_dat = olink_exposure_rsid_dat,
      outcome_dat = outcome_olink_rsid,
      action = 2
    )

    cat("Olink rsID-mapped harmonised rows:", nrow(olink_harmonised_rsid), "\n")
    print(table(olink_harmonised_rsid$mr_keep))

    cat("Olink rsID-mapped harmonisation check:\n")
    olink_harmonised_rsid %>%
      select(
        SNP,
        SNP_original,
        effect_allele.exposure,
        other_allele.exposure,
        effect_allele.outcome,
        other_allele.outcome,
        beta.exposure,
        beta.outcome,
        mr_keep,
        remove,
        palindromic,
        ambiguous
      ) %>%
      as_tibble() %>%
      print(n = Inf)

    olink_results_rsid <- run_mr_safely(
      harmonised_dat = olink_harmonised_rsid,
      label = "Olink rsID-mapped",
      harmonised_path = file.path(olink_dir, "olink_gpnmb_pd_rsid_mapped_harmonised.csv"),
    odds_path = file.path(results_dir, "olink_gpnmb_pd_rsid_mapped_odds_ratios.csv")
    )
  } else {
    message("Olink rsID-mapped MR skipped: no Parkinson's outcome rows found after rsID mapping.")
  }
} else {
  message("Olink rsID-mapped MR skipped: no Olink variants could be mapped to rsIDs.")
}

# ==============================================================================
# 7. Combined summary
# ==============================================================================

summary_results <- bind_rows(
  summarise_mr_result(
    result_object = somascan_results,
    platform = "SomaScan",
    exposure_label = "GPNMB",
    outcome_label = pd_outcome_label,
    outcome_id = pd_outcome_id
  ),
  summarise_mr_result(
    result_object = olink_results_rsid,
    platform = "Olink UKB-PPP rsID-mapped",
    exposure_label = "GPNMB",
    outcome_label = pd_outcome_label,
    outcome_id = pd_outcome_id
  )
)

if (nrow(summary_results) > 0) {
  cat("Combined MR summary:\n")
  print(summary_results)
} else {
  message("No valid MR results were generated.")
}

# ==============================================================================
# 8. Final interpretation notes
# ==============================================================================

cat("\nFinal notes:\n")
cat("- SomaScan should run if rs5850 is available in the PD outcome and harmonises cleanly.\n")
cat("- Olink chr_pos_ref_alt outcome extraction is expected to return 0 rows unless the outcome database recognises that format.\n")
cat("- Olink rsID-mapped MR depends on successful Ensembl rsID mapping and OpenGWAS outcome availability.\n")
cat("- Any single-SNP result is a Wald ratio only; do not claim robust pleiotropy or heterogeneity sensitivity analysis.\n")

# SomaScan and Olink MR analyses both suggest a positive association between genetically predicted GPNMB protein abundance and Parkinson’s disease risk. 
# The SomaScan estimate was based on a single instrument, rs5850, and therefore uses a Wald ratio. The Olink rsID-mapped analysis retained 17 harmonised 
# cis-pQTL instruments; the IVW estimate was directionally concordant with SomaScan and statistically significant. However, because the Olink variants are 
# located within the same cis locus and may be correlated, LD structure and colocalisation must be assessed before making a strong causal interpretation.