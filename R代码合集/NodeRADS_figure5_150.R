# Figure 5: ROC curves for ALNM in the 150-patient primary cohort.
# Run after installing readxl with the original Excel workbook alongside
# NodeRADS_analysis.R: Rscript NodeRADS_figure5_150.R
# Or set NODE_INPUT_XLSX, NODE_OUTPUT_DIR and NODE_ANALYSIS_R as needed.
# NODE_INPUT_CSV is an optional exact-cell export for R environments lacking
# readxl. All statistics and graphics are calculated in R.
# Raw ipsilateral/contralateral/tumor ADC are NOT candidate predictors. The
# peritumoral variable is the updated/new measurement in Excel column Z.

analysis_script <- Sys.getenv("NODE_ANALYSIS_R", unset = "NodeRADS_analysis.R")
if (!file.exists(analysis_script)) {
  stop("Analysis script not found: ", analysis_script)
}
invisible(capture.output(source(analysis_script, local = TRUE)))
stopifnot(nrow(primary) == 150L,
          sum(primary$alnm == 1) == 67L,
          sum(primary$alnm == 0) == 83L,
          setequal(names(coef(final_model))[-1],
                   c("cADC", "peri_adc_new", "size_mm", "node_final")))

out_dir <- Sys.getenv("NODE_OUTPUT_DIR", unset = "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

empirical_roc <- function(outcome, risk_score) {
  stopifnot(length(outcome) == length(risk_score),
            all(outcome %in% c(0, 1)), all(is.finite(risk_score)))
  thresholds <- c(Inf, sort(unique(risk_score), decreasing = TRUE))
  tpr <- vapply(thresholds,
                function(t) mean(risk_score[outcome == 1] >= t), numeric(1))
  fpr <- vapply(thresholds,
                function(t) mean(risk_score[outcome == 0] >= t), numeric(1))
  roc <- auc_delong(outcome, risk_score)
  auc_trapezoid <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)
  stopifnot(abs(roc["auc"] - auc_trapezoid) < 1e-12)
  list(points = data.frame(threshold = thresholds,
                           false_positive_rate = fpr,
                           sensitivity = tpr),
       auc = unname(roc["auc"]))
}

# Reverse cADC so that increasing plotted scores always mean higher ALNM risk.
# Original Node-RADS score means the adjudicated patient-level score here.
scores <- list(
  "Combined model" = predict(final_model, type = "response"),
  "Node-RADS + cADC" = predict(baseline_model, type = "response"),
  "Node-RADS" = primary$node_final,
  "cADC" = -primary$cADC,
  "Peritumoral ADC" = primary$peri_adc_new,
  "Maximum tumor diameter" = primary$size_mm
)
rocs <- lapply(scores, function(s) empirical_roc(primary$alnm, s))
auc_values <- vapply(rocs, function(x) x$auc, numeric(1))
stopifnot(identical(round(unname(auc_values), 3),
                    c(0.929, 0.888, 0.840, 0.851, 0.731, 0.755)),
          abs(auc_values["Combined model"] - full_auc["auc"]) < 1e-12,
          abs(auc_values["Node-RADS + cADC"] - base_auc["auc"]) < 1e-12)

write.csv(data.frame(Indicator = names(auc_values),
                     AUC_exact = unname(auc_values),
                     AUC_displayed = sprintf("%.3f", auc_values)),
          file.path(out_dir, "Figure5_AUC_150.csv"), row.names = FALSE)
points <- do.call(rbind, lapply(names(rocs), function(label) {
  data.frame(Indicator = label, rocs[[label]]$points)
}))
write.csv(points, file.path(out_dir, "Figure5_ROC_coordinates_150.csv"),
          row.names = FALSE)

draw_figure5 <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mar = c(4.5, 4.6, 1.8, 0.8), mgp = c(2.55, 0.65, 0),
      tcl = -0.28, las = 1, xaxs = "i", yaxs = "i",
      family = "sans", bg = "white")
  plot(NA_real_, xlim = c(-0.04, 1.04), ylim = c(0, 1),
       xlab = "", ylab = "", axes = FALSE, frame.plot = FALSE)
  abline(a = 0, b = 1, col = "#D9DDE1", lty = 2, lwd = 0.85)
  colors <- c("#FF4E3A", "#333333", "#4888D5", "#00BCA5",
              "#9F6BFA", "#FFBC32")
  widths <- c(2.0, 1.65, 1.5, 1.5, 1.5, 1.5)
  styles <- c(1, 2, 1, 1, 1, 1)
  # Draw main combined-model curve last so it stays visible at overlaps.
  for (i in rev(seq_along(rocs))) {
    path <- rocs[[i]]$points
    lines(path$false_positive_rate, path$sensitivity, type = "s",
          col = colors[i], lty = styles[i], lwd = widths[i])
  }
  ticks <- seq(0, 1, 0.2)
  axis(1, at = ticks, labels = sprintf("%.1f", ticks),
       col = "#666666", col.axis = "#222222",
       lwd = 0, lwd.ticks = 0.75, cex.axis = 0.96)
  axis(2, at = ticks, labels = sprintf("%.1f", ticks),
       col = "#666666", col.axis = "#222222",
       lwd = 0, lwd.ticks = 0.75, cex.axis = 0.96)
  box(col = "#9E9E9E", lwd = 0.75)
  mtext("1 – Specificity", side = 1, line = 2.65, cex = 1.04)
  mtext("Sensitivity", side = 2, line = 2.80, las = 0, cex = 1.04)
  mtext(sprintf("Development cohort, n = %d", nrow(primary)),
        side = 3, line = 0.38, adj = 0,
        cex = 0.85, col = "#555555")
  legend("bottomright", inset = c(0.006, 0.018), bty = "n",
         legend = sprintf("%s (AUC = %.3f)", names(auc_values), auc_values),
         col = colors, lty = styles, lwd = widths,
         cex = 0.79, seg.len = 1.9, x.intersp = 0.60,
         y.intersp = 1.02)
}

png(file.path(out_dir, "Figure_5_ROC_150.png"),
    width = 6.6, height = 6.6, units = "in", res = 360,
    pointsize = 13, type = "cairo-png")
draw_figure5()
dev.off()

tiff(file.path(out_dir, "Figure_5_ROC_150.tiff"),
     width = 6.6, height = 6.6, units = "in", res = 300,
     pointsize = 13, compression = "lzw", type = "cairo")
draw_figure5()
dev.off()

pdf(file.path(out_dir, "Figure_5_ROC_150.pdf"),
    width = 6.6, height = 6.6, pointsize = 13, useDingbats = FALSE)
draw_figure5()
dev.off()

cat("Figure 5 drawn in R for", nrow(primary), "patients; ALNM+",
    sum(primary$alnm == 1), "; ALNM-", sum(primary$alnm == 0), "\n")
print(data.frame(Indicator = names(auc_values), AUC = auc_values),
      row.names = FALSE, digits = 7)
