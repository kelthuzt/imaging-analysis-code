# Reproducible supplementary analyses for the 150-patient Node-RADS study.
# Requires NodeRADS_analysis.R and the original study workbook.
# Rscript NodeRADS_supplement_150.R
# Set NODE_INPUT_XLSX to the workbook path and NODE_OUTPUT_DIR to a results
# directory if they differ from the defaults in NodeRADS_analysis.R.
# NODE_INPUT_CSV is an optional lossless cell export used only where readxl
# is unavailable; the statistical analysis and figures below all run in R.

analysis_script <- Sys.getenv("NODE_ANALYSIS_R", "NodeRADS_analysis.R")
if (!file.exists(analysis_script)) stop("Missing analysis script: ", analysis_script)
invisible(capture.output(source(analysis_script, local = TRUE)))
stopifnot(nrow(primary) == 150L, sum(primary$alnm == 1) == 67L,
          sum(primary$alnm == 0) == 83L,
          setequal(names(coef(final_model))[-1],
                   c("node_final", "cADC", "peri_adc_new", "size_mm")))
out_dir <- Sys.getenv("NODE_OUTPUT_DIR", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# S1: all 150 participants are retained in denominators for availability.
lvi_levels <- c("LVI absent", "LVI present", "LVI unavailable")
lvi_label <- ifelse(is.na(primary$lvi), lvi_levels[3],
                    ifelse(primary$lvi == 1, lvi_levels[2], lvi_levels[1]))
subtype_levels <- c("Luminal A-like", "Luminal B-like",
                    "HER2-positive nonluminal", "Triple-negative",
                    "Subtype unavailable")
subtype_map <- c(LuminalA = subtype_levels[1],
                 LuminalB = subtype_levels[2],
                 `Her2阳性型` = subtype_levels[3], TNBC = subtype_levels[4])
subtype_label <- unname(subtype_map[primary$subtype])
subtype_label[is.na(subtype_label)] <- subtype_levels[5]
stopifnot(all(primary$subtype[!is.na(primary$subtype)] %in% names(subtype_map)))

tab_s1 <- do.call(rbind, c(
  lapply(lvi_levels, function(z) data.frame(
    Characteristic = z, Total = sum(lvi_label == z),
    Non_ALNM = sum(lvi_label == z & primary$alnm == 0),
    ALNM = sum(lvi_label == z & primary$alnm == 1))),
  lapply(subtype_levels, function(z) data.frame(
    Characteristic = z, Total = sum(subtype_label == z),
    Non_ALNM = sum(subtype_label == z & primary$alnm == 0),
    ALNM = sum(subtype_label == z & primary$alnm == 1)))
))
stopifnot(sum(tab_s1$Total[1:3]) == 150L,
          sum(tab_s1$Total[4:8]) == 150L,
          sum(!is.na(primary$lvi)) == 131L,
          sum(!is.na(primary$subtype)) == 147L)
write.csv(tab_s1, file.path(out_dir, "Supplement_Table_S1_150.csv"),
          row.names = FALSE)

# S2: LVI is postoperative. Refit the four imaging terms in the SAME complete
# subset before appending LVI; neither model reselects imaging predictors.
# This within-subset comparison cannot assess a preoperative LVI predictor.
lvi_data <- primary[!is.na(primary$lvi), , drop = FALSE]
stopifnot(nrow(lvi_data) == 131L, sum(lvi_data$alnm == 1) == 53L,
          sum(lvi_data$alnm == 0) == 78L,
          all(lvi_data$lvi %in% c(0, 1)))
imaging_131 <- glm(alnm ~ node_final + cADC + peri_adc_new + size_mm,
                   data = lvi_data, family = binomial())
extended_131 <- glm(alnm ~ node_final + cADC + peri_adc_new + size_mm + lvi,
                    data = lvi_data, family = binomial())
p_imaging <- predict(imaging_131, type = "response")
p_extended <- predict(extended_131, type = "response")
auc_imaging <- auc_delong(lvi_data$alnm, p_imaging)
auc_extended <- auc_delong(lvi_data$alnm, p_extended)
delong_lvi <- paired_auc_delong(lvi_data$alnm, p_extended, p_imaging)
lr_lvi <- anova(imaging_131, extended_131, test = "Chisq")
delta_ci <- auc_bootstrap_difference(lvi_data$alnm, p_extended,
                                     p_imaging, B = 1000L, seed = 20260929L)
scales <- c(node_final = 1, cADC = .1, peri_adc_new = .1,
            size_mm = 1, lvi = 1)
terms <- names(scales)
tab_s2 <- do.call(rbind, lapply(terms, function(term) {
  result <- wald_or(extended_131, term, scales[[term]])
  data.frame(Predictor = term, Increment = scales[[term]],
             OR = unname(result["OR"]),
             CI_lower = unname(result["CI_lower"]),
             CI_upper = unname(result["CI_upper"]),
             P = unname(result["P"]))
}))
perf_s2 <- data.frame(Model = c("Refitted four-variable imaging model",
                                "Imaging model plus postoperative LVI"),
                      N = 131L, ALNM_positive = 53L,
                      AUC = c(auc_imaging["auc"], auc_extended["auc"]),
                      CI_lower = c(auc_imaging["ci_lower"],
                                   auc_extended["ci_lower"]),
                      CI_upper = c(auc_imaging["ci_upper"],
                                   auc_extended["ci_upper"]))
comparison_s2 <- data.frame(
  AUC_difference = unname(delong_lvi["difference"]),
  Paired_DeLong_P = unname(delong_lvi["P"]),
  Bootstrap_delta_CI_lower = delta_ci[1],
  Bootstrap_delta_CI_upper = delta_ci[2],
  Likelihood_ratio_chisq = lr_lvi$Deviance[2],
  Likelihood_ratio_df = lr_lvi$Df[2],
  Likelihood_ratio_P = lr_lvi$`Pr(>Chi)`[2])
write.csv(tab_s2, file.path(out_dir, "Supplement_Table_S2_coefficients_131.csv"),
          row.names = FALSE)
write.csv(perf_s2, file.path(out_dir, "Supplement_Table_S2_performance_131.csv"),
          row.names = FALSE)
write.csv(comparison_s2, file.path(out_dir, "Supplement_S2_comparison_131.csv"),
          row.names = FALSE)

# S3: use the 150-patient imaging equation without refitting or selecting in
# each molecular subtype. Original final Node-RADS score is the comparator.
# Four paired AUC comparisons are Holm-adjusted. Small strata are exploratory.
primary$combined_probability <- predict(final_model, newdata = primary,
                                         type = "response")
primary$surrogate_subtype <- subtype_label
subtype_data <- primary[!is.na(primary$subtype), , drop = FALSE]
stopifnot(nrow(subtype_data) == 147L)
subtype_results <- do.call(rbind, lapply(subtype_levels[1:4], function(group) {
  d <- subtype_data[subtype_data$surrogate_subtype == group, , drop = FALSE]
  stopifnot(length(unique(d$alnm)) == 2L)
  im <- auc_delong(d$alnm, d$combined_probability)
  nr <- auc_delong(d$alnm, d$node_final)
  pair <- paired_auc_delong(d$alnm, d$combined_probability, d$node_final)
  data.frame(Subtype = group, N = nrow(d), Positive = sum(d$alnm == 1),
             Negative = sum(d$alnm == 0), Imaging_AUC = im["auc"],
             Imaging_CI_lower = im["ci_lower"],
             Imaging_CI_upper = im["ci_upper"],
             NodeRADS_AUC = nr["auc"],
             NodeRADS_CI_lower = nr["ci_lower"],
             NodeRADS_CI_upper = nr["ci_upper"],
             AUC_difference = pair["difference"],
             DeLong_P = pair["P"])
}))
subtype_results$Holm_P <- p.adjust(subtype_results$DeLong_P, "holm")
write.csv(subtype_results,
          file.path(out_dir, "Supplement_Table_S3_subtypes_147.csv"),
          row.names = FALSE)

roc_coords <- function(y, score) {
  thresholds <- c(Inf, sort(unique(score), decreasing = TRUE))
  data.frame(fpr = vapply(thresholds,
                          function(t) mean(score[y == 0] >= t), numeric(1)),
             tpr = vapply(thresholds,
                          function(t) mean(score[y == 1] >= t), numeric(1)))
}
rocs <- list(roc_coords(lvi_data$alnm, p_imaging),
             roc_coords(lvi_data$alnm, p_extended))
draw_s1 <- function() {
  par(mar = c(4.8, 5.0, 1.3, 1.2), mgp = c(2.8, .7, 0),
      las = 1, family = "sans", xaxs = "i", yaxs = "i")
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - Specificity",
       ylab = "Sensitivity", axes = FALSE)
  axis(1, at = seq(0, 1, .2)); axis(2, at = seq(0, 1, .2)); box()
  abline(0, 1, col = "#BBC0C4", lty = 3)
  lines(rocs[[1]]$fpr, rocs[[1]]$tpr, type = "s", lwd = 3,
        col = "#176B9A")
  lines(rocs[[2]]$fpr, rocs[[2]]$tpr, type = "s", lwd = 3,
        col = "#B45721")
  legend("bottomright", bty = "n", cex = .85, lwd = 3,
         col = c("#176B9A", "#B45721"),
         legend = c(sprintf("Imaging model (AUC %.3f)", auc_imaging["auc"]),
                    sprintf("Imaging + LVI (AUC %.3f)",
                            auc_extended["auc"])))
}

draw_s2 <- function() {
  par(mar = c(4.3, 10.2, .8, 1), mgp = c(2.5, .65, 0),
      family = "sans", las = 1)
  ys <- c(4, 3, 2, 1)
  plot(NA, xlim = c(.45, 1.05), ylim = c(.4, 4.6),
       xlab = "Area under the ROC curve", ylab = "", axes = FALSE,
       frame.plot = FALSE)
  abline(v = seq(.5, 1, .1), col = "#E6E6E6")
  axis(1, at = seq(.5, 1, .1), labels = sprintf("%.1f", seq(.5, 1, .1)))
  compact_labels <- c("Luminal A-like", "Luminal B-like",
                      "HER2-positive", "Triple-negative")
  axis(2, at = ys, labels = sprintf("%s (n=%d)",
                                   compact_labels, subtype_results$N),
       tick = FALSE, hadj = 1)
  for (i in seq_len(4)) {
    segments(subtype_results$Imaging_CI_lower[i], ys[i]+.12,
             subtype_results$Imaging_CI_upper[i], ys[i]+.12,
             col = "#176B9A", lwd = 2.5)
    segments(subtype_results$NodeRADS_CI_lower[i], ys[i]-.12,
             subtype_results$NodeRADS_CI_upper[i], ys[i]-.12,
             col = "#B45721", lwd = 2.5)
    points(subtype_results$Imaging_AUC[i], ys[i]+.12,
           pch = 16, col = "#176B9A", cex = 1.3)
    points(subtype_results$NodeRADS_AUC[i], ys[i]-.12,
           pch = 17, col = "#B45721", cex = 1.3)
  }
  legend("bottomright", bty = "n", cex = .9, pch = c(16, 17),
         col = c("#176B9A", "#B45721"),
         legend = c("Combined imaging model", "Node-RADS"))
}

for (figure in c("S1", "S2")) {
  dest <- file.path(out_dir, paste0("Figure_", figure, "_150.png"))
  png(dest, width = 7.2, height = 5.4, units = "in", res = 300,
      pointsize = 14, type = "cairo-png")
  if (figure == "S1") draw_s1() else draw_s2()
  dev.off()
}

cat("Supplementary cohort:", nrow(primary), "total,", nrow(lvi_data),
    "LVI observed,", nrow(subtype_data), "subtype observed\n")
print(tab_s1, row.names = FALSE)
print(tab_s2, row.names = FALSE)
print(perf_s2, row.names = FALSE)
print(comparison_s2, row.names = FALSE)
print(subtype_results, row.names = FALSE)
