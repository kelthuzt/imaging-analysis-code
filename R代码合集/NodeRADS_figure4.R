# Figure 4: Node-RADS score pairs, cross-threshold adjudication and nodal ADC.
# Run this in the same directory as NodeRADS_analysis.R after installing readxl:
#   Rscript NodeRADS_figure4.R
# Set NODE_INPUT_XLSX and NODE_OUTPUT_DIR if the data/output live elsewhere.
# NODE_INPUT_CSV is an optional deidentified exact-cell fallback for R runtimes
# without readxl. All plotted counts, ADC values and P labels come from data.

analysis_script <- Sys.getenv("NODE_ANALYSIS_R", unset = "NodeRADS_analysis.R")
if (!file.exists(analysis_script)) stop("Analysis script not found: ", analysis_script)
invisible(capture.output(source(analysis_script, local = TRUE)))
stopifnot(nrow(primary) == 150L, sum(primary$alnm == 1L) == 67L)

output_dir <- Sys.getenv("NODE_OUTPUT_DIR", unset = "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

score_levels <- 1:5
score_pairs <- table(
  Senior = factor(primary$node_senior, levels = score_levels),
  Junior = factor(primary$node_junior, levels = score_levels)
)
senior_pred <- primary$node_senior >= 3
junior_pred <- primary$node_junior >= 3
cross_threshold <- senior_pred != junior_pred
positive <- primary$alnm == 1
correct_cross <- c(
  Senior = sum(senior_pred[cross_threshold] == positive[cross_threshold]),
  Junior = sum(junior_pred[cross_threshold] == positive[cross_threshold]),
  Adjudicated = sum((primary$node_final[cross_threshold] >= 3) ==
                     positive[cross_threshold])
)
cross_n <- sum(cross_threshold)
pair_p <- binom.test(correct_cross["Senior"],
                     correct_cross["Senior"] + correct_cross["Junior"],
                     p = 0.5)$p.value

both_missed <- positive & !senior_pred & !junior_pred
both_detected <- positive & senior_pred & junior_pred
adc_missed <- primary$ipsi_adc[both_missed]
adc_detected <- primary$ipsi_adc[both_detected]
adc_holm_p <- missed_results$Holm_P[missed_results$Measure == "ipsi_adc"]
stopifnot(sum(score_pairs) == nrow(primary), cross_n == 14L,
          sum(both_missed) == 18L, sum(both_detected) == 44L,
          length(adc_holm_p) == 1L)

draw_figure4 <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  layout(matrix(1:3, nrow = 1), widths = c(1.06, 1, 1.08))
  par(bg = "white", family = "sans", cex.axis = 0.98, cex.lab = 1.02,
      mgp = c(1.95, 0.50, 0), tck = -0.018, las = 1)

  # A. Individual five-level score pairs. Red boxes mark cells for which
  # classification changes across the prespecified Node-RADS threshold of 3.
  par(mar = c(3.6, 3.7, 2.7, 0.30))
  plot(NA_real_, xlim = c(0.5, 5.5), ylim = c(5.5, 0.5),
       xlab = "Junior reader score", ylab = "Senior reader score",
       axes = FALSE, asp = 1, xaxs = "i", yaxs = "i")
  blue_scale <- colorRampPalette(c("#FFFFFF", "#D8E7F2", "#265678"))(
    max(score_pairs) + 1L)
  for (i in score_levels) for (j in score_levels) {
    value <- score_pairs[i, j]
    rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5,
         col = blue_scale[value + 1L], border = "#B9C5CD", lwd = 0.7)
    text(j, i, labels = value, cex = 1.02,
         col = if (value >= 35L) "white" else "#28343D")
  }
  rect(2.5, 0.5, 5.5, 2.5, border = "#C74240", lwd = 1.5)
  rect(0.5, 2.5, 2.5, 5.5, border = "#C74240", lwd = 1.5)
  axis(1, at = score_levels, labels = score_levels)
  axis(2, at = score_levels, labels = score_levels)
  box(col = "#677782", lwd = 0.8)
  title(main = "A  Reader score pairs", adj = 0, line = 0.40,
        cex.main = 1.02, font.main = 2)

  # B. Paired classification of the same cross-threshold patients.
  par(mar = c(3.6, 3.3, 2.7, 0.25))
  ylim_top <- max(14, max(correct_cross) + 2)
  bar_x <- barplot(correct_cross, ylim = c(0, ylim_top),
                   names.arg = rep("", 3),
                   col = c("#25577D", "#58A640", "#ED822C"),
                   border = NA, width = 0.72, space = 0.42,
                   ylab = "Correct classifications", axes = FALSE)
  axis(2, at = seq(0, ylim_top, by = 2))
  axis(1, at = bar_x, labels = c("Senior\nreader", "Junior\nreader",
                               "Third\nscore"), tick = FALSE, line = 1)
  box(bty = "l", col = "#677782", lwd = 0.8)
  text(bar_x, correct_cross + 0.45,
       labels = paste0(correct_cross, "/", cross_n), cex = 0.99)
  text(mean(range(bar_x)), ylim_top - 0.32,
       sprintf("Senior vs junior: P = %.3f", pair_p), cex = 0.86)
  title(main = "B  Cross-threshold cases",
        adj = 0, line = 0.40, cex.main = 1.02, font.main = 2)

  # C. Match the reference: filled interquartile-range rectangles with
  # individual patient dots, without whiskers or a median line. The ADC unit
  # and group denominators are specified in the figure caption.
  # Reinitialize the jitter seed for identical dots in PDF/PNG/TIFF exports.
  par(mar = c(3.6, 3.5, 2.7, 0.25))
  plot(NA_real_, xlim = c(0.48, 2.52), ylim = c(0.50, 2.60),
       xlab = "", ylab = "Axillary node ADC", axes = FALSE,
       xaxs = "i", yaxs = "i")
  q_missed <- quantile(adc_missed, c(0.25, 0.75), names = FALSE)
  q_detected <- quantile(adc_detected, c(0.25, 0.75), names = FALSE)
  rect(1 - 0.30, q_missed[1], 1 + 0.30, q_missed[2],
       col = "#CF5757", border = NA)
  rect(2 - 0.30, q_detected[1], 2 + 0.30, q_detected[2],
       col = "#477EC0", border = NA)
  axis(2, at = seq(0.5, 2.5, by = 0.5))
  axis(1, at = 1:2, labels = c("Both missed", "Both detected"),
       tick = FALSE, line = 1)
  box(bty = "l", col = "#677782", lwd = 0.8)
  set.seed(20260928L)
  points(jitter(rep(1, length(adc_missed)), amount = 0.13), adc_missed,
         pch = 16, cex = 0.35, col = adjustcolor("#26343D", alpha.f = 0.75))
  points(jitter(rep(2, length(adc_detected)), amount = 0.13), adc_detected,
         pch = 16, cex = 0.35, col = adjustcolor("#26343D", alpha.f = 0.75))
  p_label <- if (adc_holm_p < 0.001) "Holm-adjusted P < 0.001" else
    sprintf("Holm-adjusted P = %.3f", adc_holm_p)
  text(1.5, 2.47, p_label, cex = 0.86)
  title(main = "C  ALNM-positive cases", adj = 0, line = 0.40,
        cex.main = 1.02, font.main = 2)
}

png(file.path(output_dir, "Figure_4_NodeRADS_150.png"),
    width = 10.5, height = 4.1, units = "in", res = 300,
    pointsize = 18, type = "cairo-png")
draw_figure4()
dev.off()

tiff(file.path(output_dir, "Figure_4_NodeRADS_150.tiff"),
     width = 10.5, height = 4.1, units = "in", res = 300,
     pointsize = 18, compression = "lzw", type = "cairo")
draw_figure4()
dev.off()

pdf(file.path(output_dir, "Figure_4_NodeRADS_150.pdf"),
    width = 10.5, height = 4.1, pointsize = 18, useDingbats = FALSE)
draw_figure4()
dev.off()

cat("Figure 4 saved; cohort n =", nrow(primary),
    "; cells (1,1) and (3,3) =", score_pairs[1, 1], score_pairs[3, 3],
    "; cross-threshold =", cross_n,
    "; jointly missed/detected =", length(adc_missed),
    length(adc_detected), "\n")
