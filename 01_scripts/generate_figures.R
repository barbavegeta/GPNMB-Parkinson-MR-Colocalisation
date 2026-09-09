#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# Dissertation figures 2 to 5.
#
# Base graphics only - no package dependencies.
#
# Usage:
#   Rscript generate_figures.R              # -> 02_results/figures/
#   Rscript generate_figures.R path/to/out
# ---------------------------------------------------------------------

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit)) return(normalizePath(sub("^--file=", "", hit[1])))
  if (!is.null(sys.frames()[[1]]$ofile)) return(normalizePath(sys.frames()[[1]]$ofile))
  NA_character_
}

sp <- script_path()
HERE <- if (!is.na(sp)) dirname(sp) else getwd()
ROOT <- if (basename(HERE) == "01_scripts") dirname(HERE) else getwd()
if (!dir.exists(file.path(ROOT, "02_results"))) {
  stop(paste("Could not find 02_results/ under", ROOT,
             "- run this from the package root or from 01_scripts/."))
}

cli <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(cli)) cli[1] else file.path(ROOT, "02_results", "figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

read_required <- function(relative_path, glob_fallback = NULL) {
  path <- file.path(ROOT, relative_path)
  if (!file.exists(path) && !is.null(glob_fallback)) {
    hits <- sort(Sys.glob(file.path(ROOT, dirname(relative_path), glob_fallback)))
    if (length(hits)) path <- hits[1]
  }
  if (!file.exists(path)) stop(sprintf("Required input not found: %s", path))
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

# ---------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------
SHOW_TITLES <- FALSE     # captions carry the title in the dissertation
DPI         <- 300

TEAL  <- "#005B5E"; NAVY  <- "#001A79"; STEEL <- "#2F6E8A"
LBLUE <- "#7BB9CB"; GREY  <- "#5A6470"; ORANGE <- "#C2410C"
INK   <- "#1B2129"; MUTED <- "#3D4650"
GRID  <- "#E7EBEE"; POINT <- "#9FB6C4"; AXIS <- "#9AA4AC"

pt   <- function(points) points * 96 / 72
mcex <- function(d_pt) (d_pt - 0.96) / 4.32

open_png <- function(file, w, h, ps = 9) {
  png(file.path(OUT, file), width = w, height = h, units = "in",
      res = DPI, type = "cairo", bg = "white", pointsize = ps)
}
done <- function(file) { dev.off(); message(sprintf("Wrote: %s", file.path(OUT, file))) }

ax <- function(side, at, labels = TRUE, cex = 1, tcl = -0.3, ...) {
  axis(side, at = at, labels = labels, tcl = tcl, lwd = 0,
       lwd.ticks = pt(0.8), col.ticks = AXIS, col.axis = MUTED,
       cex.axis = cex, ...)
}
spine <- function(side, from, to) {
  u <- par("usr")
  if (side == 1) segments(from, u[3], to, u[3], col = AXIS, lwd = pt(0.9), xpd = NA)
  if (side == 2) segments(u[1], from, u[1], to, col = AXIS, lwd = pt(0.9), xpd = NA)
}
maybe_title <- function(txt) if (SHOW_TITLES) mtext(txt, side = 3, line = 0.6, cex = 1.15, col = INK)


# =====================================================================
# Figure 2 - GPNMB RNA expression across brain regions
# =====================================================================
# Human Protein Atlas consensus values as reported in Appendix I.
# The choroid plexus is a CSF-producing epithelial and vascular structure
# in the ventricles, not neural tissue, so it is shaded separately and
# excluded from the neural ranking given in the text.
hpa <- data.frame(
  tissue = c("choroid plexus", "white matter", "spinal cord", "pons", "thalamus",
             "medulla oblongata", "amygdala", "cerebral cortex", "basal ganglia",
             "midbrain", "hippocampal formation", "pituitary gland",
             "hypothalamus", "cerebellum", "retina"),
  nTPM = c(146.4, 109.6, 49.4, 35.8, 34.6, 32.2, 29.4, 24.4, 16.5,
           14.5, 14.3, 12.7, 8.9, 3.2, 2.6),
  neural = c(FALSE, rep(TRUE, 14)),
  stringsAsFactors = FALSE
)
hpa <- hpa[order(hpa$nTPM), ]

open_png("Figure_2_GPNMB_brain_RNA_expression.png", 6.9, 5.0)
par(mai = c(0.52, 1.78, 0.30, 0.22), xaxs = "i", yaxs = "i", las = 1)

n <- nrow(hpa); yc <- seq_len(n) - 1
plot(NA, xlim = c(0, 158), ylim = c(-0.8, n - 1 + 0.8),
     axes = FALSE, xlab = "", ylab = "")
abline(v = seq(0, 140, by = 20), col = GRID, lwd = pt(0.8))
rect(0, yc - 0.4, hpa$nTPM, yc + 0.4,
     col = ifelse(hpa$neural, TEAL, LBLUE), border = NA)
text(hpa$nTPM + 2.5, yc, format(hpa$nTPM, nsmall = 1),
     adj = c(0, 0.5), col = MUTED, cex = 0.86)

ax(1, seq(0, 140, by = 20), cex = 0.95)
ax(2, yc, hpa$tissue, cex = 0.95)
spine(1, 0, 158); spine(2, -0.8, n - 1 + 0.8)
mtext("Consensus RNA expression (nTPM)", side = 1, line = 1.9, col = INK, cex = 1.0)
mtext("Brain region", side = 2, line = 10.3, las = 0, col = INK, cex = 1.0)
legend("bottomright", inset = c(0.015, 0.02),
       legend = c("Neural tissue", "Non-neural (CSF-producing epithelium)"),
       fill = c(TEAL, LBLUE), border = NA, bty = "n",
       cex = 0.86, text.col = MUTED, y.intersp = 1.15)
maybe_title("GPNMB RNA expression across brain regions")
done("Figure_2_GPNMB_brain_RNA_expression.png")


# =====================================================================
# Figure 3 - Mendelian randomisation estimates
# =====================================================================
soma <- read_required(
  "02_results/current_results/somascan/somascan_gpnmb_pd_odds_ratios.csv")[1, ]
lead <- read_required(
  "02_results/latest_olink/gpnmb_corrected_mr_outputs/01_lead_variant_wald_mr.csv")[1, ]
fixed <- read_required(
  "02_results/latest_olink/gpnmb_ld_clumped_mr_outputs/05_all_ld_clumped_mr_results.csv",
  glob_fallback = "*_all_ld_clumped_mr_results.csv")
rownames(fixed) <- fixed$setting
re_ivw <- read_required(
  "02_results/latest_olink/gpnmb_ld_clumped_mr_outputs/06_multiplicative_random_effects_ivw_results.csv",
  glob_fallback = "*_multiplicative_random_effects_ivw_results.csv")
rownames(re_ivw) <- re_ivw$setting
st <- "strict_r2_0p001_10Mb"; sd_ <- "standard_r2_0p01_10Mb"

# Colour encodes the type of analysis, not merely the row.
forest <- data.frame(
  label = c("SomaScan, Wald ratio (1 SNP)",
            "Olink, lead variant rs75801644",
            "Olink, strict clumping, fixed IVW (5 SNPs)",
            "Olink, strict clumping, random effects",
            "Olink, standard clumping, fixed IVW (12 SNPs)",
            "Olink, standard clumping, random effects"),
  est = c(soma$or, lead$mr_OR, fixed[st, "mr_OR"], re_ivw[st, "multiplicative_RE_OR"],
          fixed[sd_, "mr_OR"], re_ivw[sd_, "multiplicative_RE_OR"]),
  lo  = c(soma$or_lci95, lead$mr_OR_lci, fixed[st, "mr_OR_lci"],
          re_ivw[st, "multiplicative_RE_OR_lci"], fixed[sd_, "mr_OR_lci"],
          re_ivw[sd_, "multiplicative_RE_OR_lci"]),
  hi  = c(soma$or_uci95, lead$mr_OR_uci, fixed[st, "mr_OR_uci"],
          re_ivw[st, "multiplicative_RE_OR_uci"], fixed[sd_, "mr_OR_uci"],
          re_ivw[sd_, "multiplicative_RE_OR_uci"]),
  col = c(NAVY, GREY, TEAL, STEEL, TEAL, STEEL),
  stringsAsFactors = FALSE
)

# p-values are drawn with plotmath rather than Unicode superscripts, which
# are missing from some system fonts and render as an empty box.
plabs <- list(bquote("2.6" %*% 10^-7), "0.165", "0.057", "0.473", "0.017", "0.221")

open_png("Figure_3_GPNMB_MR_forest_plot.png", 7.4, 4.3)
par(mai = c(0.56, 2.52, 0.34, 0.18), xaxs = "i", yaxs = "i", las = 1)

k <- nrow(forest); yy <- rev(seq_len(k)) - 1
plot(NA, xlim = c(0.80, 1.66), ylim = c(-0.75, k - 0.25),
     axes = FALSE, xlab = "", ylab = "")
xt <- seq(0.8, 1.5, by = 0.1)
abline(v = xt, col = GRID, lwd = pt(0.8))
abline(v = 1.0, col = INK, lty = 2, lwd = pt(1.1))

cap <- 0.12
for (i in seq_len(k)) with(forest[i, ], {
  lab <- plabs[[i]]
  segments(lo, yy[i], hi, yy[i], col = col, lwd = pt(2.0), lend = "butt")
  segments(lo, yy[i] - cap, lo, yy[i] + cap, col = col, lwd = pt(1.6))
  segments(hi, yy[i] - cap, hi, yy[i] + cap, col = col, lwd = pt(1.6))
  points(est, yy[i], pch = 21, bg = col, col = "white", cex = mcex(8), lwd = pt(0.9))
  text(1.645, yy[i], lab, adj = c(1, 0.5), col = INK, cex = 0.88, xpd = NA)
})
text(1.645, k - 0.45, expression(italic(p)), adj = c(1, 0.5),
     col = MUTED, cex = 0.9, xpd = NA)

ax(1, xt, format(xt, nsmall = 1), cex = 0.92)
ax(2, yy, forest$label, cex = 0.88)
spine(1, 0.80, 1.50)
mtext("Odds ratio for Parkinson's disease (95% CI)", side = 1, line = 1.8,
      at = 1.15, col = INK, cex = 1.0)
maybe_title("Mendelian randomisation estimates for genetically predicted GPNMB")
done("Figure_3_GPNMB_MR_forest_plot.png")


# =====================================================================
# Figure 4 - Regional association patterns
# =====================================================================
pqtl <- read_required(
  "02_results/current_results/regional_inputs/gpnmb_pqtl_region_summary_peer_review.csv")
pdg <- read_required(
  "02_results/current_results/regional_inputs/pd_gwas_gpnmb_region_summary_peer_review.csv")
for (nm in c("pqtl", "pdg")) {
  d <- get(nm); d$mb <- d$pos / 1e6
  d$lp <- -log10(pmax(d$p, .Machine$double.xmin)); assign(nm, d)
}

# Leader-line offsets: x nudge (Mb), y nudge (-log10 p), text adj.
hl_pqtl <- list(rs5850 = c(0.09, 1.5, 0), rs4140959 = c(-0.11, 1.6, 1),
                rs6461688 = c(-0.13, 2.1, 1))
hl_pd   <- list(rs858295 = c(0.08, 1.7, 0), rs199347 = c(0.19, 0.7, 0),
                rs6461688 = c(-0.14, 1.7, 1))

XLIM4 <- range(c(pqtl$mb, pdg$mb)) + c(-0.022, 0.022)   # so the 22.8 tick is drawn

panel <- function(d, label, hl, colour, xaxis) {
  ymax <- max(d$lp)
  plot(NA, xlim = XLIM4, ylim = c(-0.6, ymax * 1.34),
       axes = FALSE, xlab = "", ylab = "")
  yt <- pretty(c(0, ymax)); yt <- yt[yt <= ymax * 1.34]
  abline(h = yt, col = GRID, lwd = pt(0.8))
  points(d$mb, d$lp, pch = 19, cex = mcex(3.6),
         col = adjustcolor(POINT, alpha.f = 0.75))
  for (snp in names(hl)) {
    r <- d[d$snp == snp, ]; if (!nrow(r)) next
    o <- hl[[snp]]; px <- r$mb[1]; py <- r$lp[1]
    segments(px, py, px + o[1], py + o[2], col = colour, lwd = pt(0.8))
    points(px, py, pch = 21, bg = colour, col = "white", cex = mcex(7), lwd = pt(0.9))
    text(px + o[1], py + o[2], snp, adj = c(o[3], -0.3), col = colour, cex = 0.84)
  }
  ax(2, yt, cex = 0.9); spine(2, min(yt), max(yt))
  mtext(label, side = 3, line = 0.25, adj = 0, col = colour, cex = 0.95)
  mtext(expression(-log[10](italic(p))), side = 2, line = 2.0, las = 0,
        col = INK, cex = 0.95)
  if (xaxis) {
    ax(1, seq(22.8, 23.8, by = 0.2), cex = 0.9)
    spine(1, XLIM4[1], XLIM4[2])
    mtext("Chromosome 7 position (Mb, GRCh37)", side = 1, line = 1.8,
          col = INK, cex = 1.0)
  }
}

open_png("Figure_4_GPNMB_regional_association.png", 6.9, 5.4)
par(mfrow = c(2, 1), xaxs = "i", yaxs = "i", las = 1)
par(mai = c(0.22, 0.66, 0.28, 0.16))
panel(pqtl, "GPNMB protein quantitative trait locus (SomaScan, 165 variants)",
      hl_pqtl, TEAL, FALSE)
par(mai = c(0.54, 0.66, 0.28, 0.16))
panel(pdg, "Parkinson's disease GWAS (ieu-b-7, 6,368 variants)",
      hl_pd, NAVY, TRUE)
maybe_title("Regional association patterns across the GPNMB locus")
done("Figure_4_GPNMB_regional_association.png")


# =====================================================================
# Figure 5 - Colocalisation prior sensitivity
# =====================================================================
prior <- read_required("02_results/current_results/gpnmb_coloc_prior_sensitivity.csv")
prior <- prior[order(prior$p12), ]
lx <- log10(prior$p12)

open_png("Figure_5_GPNMB_coloc_prior_sensitivity.png", 6.9, 4.3)
par(mai = c(0.56, 0.66, 0.30, 0.20), xaxs = "i", yaxs = "i", las = 1)
plot(NA, xlim = c(-7.22, -3.88), ylim = c(-0.04, 1.10),
     axes = FALSE, xlab = "", ylab = "")

hg <- seq(0, 1, by = 0.2)
abline(h = hg, col = GRID, lwd = pt(0.8))
abline(v = -7:-4, col = GRID, lwd = pt(0.8))

lines(lx, prior[["PP.H4"]], col = TEAL, lwd = pt(2.0))
lines(lx, prior[["PP.H3"]], col = ORANGE, lwd = pt(2.0))
points(lx, prior[["PP.H4"]], pch = 21, bg = TEAL, col = "white",
       cex = mcex(7), lwd = pt(0.9))
points(lx, prior[["PP.H3"]], pch = 21, bg = ORANGE, col = "white",
       cex = mcex(7), lwd = pt(0.9))

# The two priors quoted in the text, plus the crossover.
cross <- approx(prior[["PP.H4"]] - prior[["PP.H3"]], lx, xout = 0)$y
crossy <- approx(lx, prior[["PP.H4"]], xout = cross)$y
segments(cross, -0.04, cross, crossy, col = MUTED, lty = 3, lwd = pt(0.9))
text(cross + 0.05, 0.045, "PP.H4 = PP.H3", adj = c(0, 0.5), col = MUTED, cex = 0.8)

points(log10(1e-5), 0.9440, pch = 21, bg = "white", col = TEAL, cex = mcex(10), lwd = pt(1.4))
text(log10(1e-5) - 0.06, 0.9440 - 0.10, "0.944", adj = c(1, 0.5), col = TEAL, cex = 0.86)
points(log10(1e-7), 0.1442, pch = 21, bg = "white", col = TEAL, cex = mcex(10), lwd = pt(1.4))
text(log10(1e-7) + 0.07, 0.1442 - 0.075, "0.144", adj = c(0, 0.5), col = TEAL, cex = 0.86)

ax(1, -7:-4, parse(text = paste0("10^", -7:-4)), cex = 0.92)
minor <- as.vector(outer(log10(2:9), -7:-4, "+"))
axis(1, at = minor[minor > -7.22 & minor < -3.88], labels = FALSE, tcl = -0.16,
     lwd = 0, lwd.ticks = pt(0.8), col.ticks = AXIS)
ax(2, hg, format(hg, nsmall = 1), cex = 0.92)
spine(1, -7.22, -3.88); spine(2, 0, 1)

mtext(expression("Shared-causal-variant prior, " * italic(p)[12]), side = 1,
      line = 1.9, col = INK, cex = 1.0)
mtext("Posterior probability", side = 2, line = 2.1, las = 0, col = INK, cex = 1.0)
legend(-7.22, 1.10, legend = c("PP.H4, shared causal variant",
                               "PP.H3, distinct causal variants"),
       col = c(TEAL, ORANGE), lwd = pt(2.0), pch = 19, pt.cex = mcex(7),
       bty = "n", cex = 0.84, text.col = MUTED, seg.len = 1.4,
       y.intersp = 1.2, horiz = FALSE)
maybe_title("Sensitivity of PP.H3 and PP.H4 to the shared-causal-variant prior")
done("Figure_5_GPNMB_coloc_prior_sensitivity.png")

message(sprintf("\nFinished. Dissertation figures saved in: %s", OUT))

