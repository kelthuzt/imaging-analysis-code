# Figure 6: internal validation, decision curves, and four-variable nomogram.
# Run from the directory holding NodeRADS_analysis.R and original XLSX:
#   install.packages("readxl")                  # only if necessary
#   Sys.setenv(NODE_INPUT_XLSX="your_data.xlsx", NODE_OUTPUT_DIR="results")
#   source("NodeRADS_figure6_150.R")
# Optional NODE_INPUT_CSV uses the exact deidentified spreadsheet cell export.
# All statistical calculations and figures are made in R; no patient rows
# are exported. R's base graphics, stats and grDevices packages suffice.

analysis_script <- Sys.getenv("NODE_ANALYSIS_R", unset = "NodeRADS_analysis.R")
if (!file.exists(analysis_script)) stop("Analysis script not found: ", analysis_script)
invisible(capture.output(source(analysis_script, local = TRUE)))
stopifnot(nrow(primary) == 150L, sum(primary$alnm == 1) == 67L,
          sum(primary$alnm == 0) == 83L)
out_dir <- Sys.getenv("NODE_OUTPUT_DIR", unset = "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Recreate the manuscript's ENTIRE model-development pipeline inside each
# bootstrap sample: all six preoperative candidate predictors, univariable
# P < .10 screening, and backward P >= .05 removal. Continuous terms use
# Wald P; TIC (three-level factor) uses a single group likelihood-ratio P.
# TIC is not estimable as a three-level factor if a bootstrap sample misses
# the rare third level; that sample sets its screening P to 1 rather than
# constructing a non-predictable model for the original patients.
continuous_candidates <- c("cADC", "rADC", "peri_adc_new", "size_mm",
                           "node_final")
model_pipeline <- function(dat) {
  p <- vapply(continuous_candidates, function(v) {
    fit <- suppressWarnings(glm(reformulate(v, response = "alnm"),
                                data = dat, family = binomial()))
    value <- coef(summary(fit))[v, "Pr(>|z|)"]
    if (is.finite(value)) value else 1
  }, numeric(1))
  p_tic <- if (length(unique(dat$tic)) < 3L) 1 else {
    null <- glm(alnm ~ 1, data = dat, family = binomial())
    test <- suppressWarnings(glm(alnm ~ factor(tic), data = dat,
                                 family = binomial()))
    value <- anova(null, test, test = "Chisq")$`Pr(>Chi)`[2]
    if (is.finite(value)) value else 1
  }
  p <- c(p, "factor(tic)" = p_tic)
  active <- names(p)[p < 0.10]
  if (!length(active)) {
    return(list(fit = glm(alnm ~ 1, data = dat, family = binomial()),
                terms = "(Intercept)"))
  }
  repeat {
    f <- as.formula(paste("alnm ~", paste(active, collapse = " + ")))
    fit <- suppressWarnings(glm(f, data = dat, family = binomial()))
    term_p <- vapply(active, function(v) {
      value <- if (v == "factor(tic)") {
        suppressWarnings(drop1(fit, test = "Chisq")[
          "factor(tic)", "Pr(>Chi)"])
      } else {
        coef(summary(fit))[v, "Pr(>|z|)"]
      }
      if (is.finite(value)) value else 1
    }, numeric(1))
    if (all(term_p < 0.05)) return(list(fit = fit, terms = active))
    active <- setdiff(active, names(which.max(term_p)))
    if (!length(active)) {
      return(list(fit = glm(alnm ~ 1, data = dat, family = binomial()),
                  terms = "(Intercept)"))
    }
  }
}

original_pipeline <- model_pipeline(primary)
stopifnot(setequal(original_pipeline$terms,
                   c("cADC", "peri_adc_new", "size_mm", "node_final")),
          all(abs(coef(original_pipeline$fit)[names(coef(final_model))] -
                    coef(final_model)) < 1e-10))

# Apparent performance (same patients used to build and evaluate the model).
apparent_probs <- predict(final_model, type = "response")
apparent_auc <- unname(auc_delong(primary$alnm, apparent_probs)["auc"])
apparent_lp <- predict(final_model, type = "link")
apparent_cal_intercept <- unname(coef(glm(
  primary$alnm ~ 1 + offset(apparent_lp), family = binomial()))[1])
apparent_cal_slope <- unname(coef(glm(
  primary$alnm ~ apparent_lp, family = binomial()))[2])
stopifnot(abs(apparent_cal_intercept) < 1e-8,
          abs(apparent_cal_slope - 1) < 1e-8)

# In-sample bootstrap performance - performance when each bootstrap model is
# applied back to the original 150 patients = optimism for that replicate.
# Logistic recalibration in each bootstrap training dataset gives intercept 0
# and slope 1 (MLE property); testing on the original data gives the test
# calibration intercept and slope. The means are the optimism-corrected
# intercept and slope of the development procedure.
B <- 1000L
bootstrap_seed <- 20260928L
set.seed(bootstrap_seed)
boot_rows <- vector("list", B)
for (b in seq_len(B)) {
  index <- sample.int(nrow(primary), replace = TRUE)
  sample_data <- primary[index, , drop = FALSE]
  selected <- model_pipeline(sample_data)
  lp_train <- suppressWarnings(predict(selected$fit,
                                       newdata = sample_data, type = "link"))
  lp_test <- suppressWarnings(predict(selected$fit,
                                      newdata = primary, type = "link"))
  if (!all(is.finite(c(lp_train, lp_test)))) {
    stop("Bootstrap iteration ", b, " has non-finite predictions")
  }
  if (sd(lp_test) <= 0) {
    stop("Bootstrap iteration ", b, " has a constant test predictor")
  }
  auc_train <- unname(auc_delong(sample_data$alnm, lp_train)["auc"])
  auc_test <- unname(auc_delong(primary$alnm, lp_test)["auc"])
  test_intercept <- unname(coef(glm(
    primary$alnm ~ 1 + offset(lp_test), family = binomial()))[1])
  test_slope <- unname(coef(glm(
    primary$alnm ~ lp_test, family = binomial()))[2])
  if (!all(is.finite(c(auc_train, auc_test, test_intercept, test_slope)))) {
    stop("Bootstrap iteration ", b, " has non-finite performance")
  }
  boot_rows[[b]] <- data.frame(
    Replicate = b, Distinct_patients = length(unique(index)),
    Selected_predictors = paste(selected$terms, collapse = ";"),
    Bootstrap_auc = auc_train, Original_data_auc = auc_test,
    Auc_optimism = auc_train - auc_test,
    Original_data_calibration_intercept = test_intercept,
    Original_data_calibration_slope = test_slope)
  if (b %% 100L == 0L) cat("Internal validation:", b, "of", B, "resamples\n")
}
boot <- do.call(rbind, boot_rows)
stopifnot(nrow(boot) == 1000L, all(complete.cases(boot)))
corrected_auc <- apparent_auc - mean(boot$Auc_optimism)
corrected_intercept <- apparent_cal_intercept -
  mean(0 - boot$Original_data_calibration_intercept)
corrected_slope <- apparent_cal_slope -
  mean(1 - boot$Original_data_calibration_slope)
validation_summary <- data.frame(
  Metric = c("Primary cohort n", "ALNM positive events", "Bootstrap resamples",
             "Random seed", "Apparent AUC", "Mean AUC optimism",
             "Optimism-corrected AUC", "Apparent calibration intercept",
             "Apparent calibration slope", "Corrected calibration intercept",
             "Corrected calibration slope"),
  Value = c(150, 67, B, bootstrap_seed, apparent_auc,
            mean(boot$Auc_optimism), corrected_auc,
            apparent_cal_intercept, apparent_cal_slope,
            corrected_intercept, corrected_slope))
write.csv(boot, file.path(out_dir, "06_bootstrap_1000_replicates.csv"),
          row.names = FALSE)
write.csv(validation_summary,
          file.path(out_dir, "06_validation_summary_150.csv"),
          row.names = FALSE)
selected_count <- as.data.frame(table(boot$Selected_predictors))
names(selected_count) <- c("Selected_predictors", "Replicates")
write.csv(selected_count,
          file.path(out_dir, "06_bootstrap_selection_frequency.csv"),
          row.names = FALSE)

# The corrected calibration line in panel A is logistic recalibration based
# on the bootstrap-corrected intercept/slope. Observed points are simple
# apparent decile event rates (15 patients per bin), NOT a claim that a smooth
# observed calibration curve was itself corrected by the bootstrap.
risk_decile <- cut(rank(apparent_probs, ties.method = "first"),
                   breaks = seq(0, 150, by = 15), include.lowest = TRUE)
calibration_bins <- do.call(rbind, lapply(split(seq_len(150), risk_decile),
                                    function(i) data.frame(
  N = length(i), Mean_predicted_risk = mean(apparent_probs[i]),
  Observed_ALNM_rate = mean(primary$alnm[i]))))
calibration_bins$Decile <- seq_len(nrow(calibration_bins))
stopifnot(nrow(calibration_bins) == 10L,
          all(calibration_bins$N == 15L))
write.csv(calibration_bins,
          file.path(out_dir, "06_calibration_deciles_150.csv"),
          row.names = FALSE)

# Decision curves based on development-cohort predicted probabilities.
node_model <- glm(alnm ~ node_final, data = primary, family = binomial())
node_probs <- predict(node_model, type = "response")
thresholds <- seq(0.01, 0.95, by = 0.01)
net_benefit <- function(y, p, t) {
  mean(y == 1 & p >= t) - mean(y == 0 & p >= t) * t / (1 - t)
}
dca <- data.frame(
  Threshold_probability = thresholds,
  Combined_model = vapply(thresholds,
    function(t) net_benefit(primary$alnm, apparent_probs, t), numeric(1)),
  Node_RADS_only = vapply(thresholds,
    function(t) net_benefit(primary$alnm, node_probs, t), numeric(1)),
  Treat_all = mean(primary$alnm == 1) -
    mean(primary$alnm == 0) * thresholds / (1 - thresholds),
  Treat_none = 0)
dca$Combined_highest <- dca$Combined_model > pmax(
  dca$Node_RADS_only, dca$Treat_all, dca$Treat_none)
range_index <- which(dca$Threshold_probability >= .10 - 1e-9 &
                     dca$Threshold_probability <= .80 + 1e-9)
stopifnot(all(dca$Combined_highest[range_index]))
write.csv(dca, file.path(out_dir, "06_decision_curve_150.csv"),
          row.names = FALSE)

# Algebraically exact nomogram: for each predictor use the least risky end of
# its observed range as 0 points, and set the largest possible contribution
# from one predictor to 100 points. Total points map back to the SAME fitted
# (unshrunk) logistic model, without rounding the model coefficients.
nom_vars <- c("node_final", "size_mm", "cADC", "peri_adc_new")
betas <- coef(final_model)[nom_vars]
ranges <- lapply(nom_vars, function(v) range(primary[[v]]))
names(ranges) <- nom_vars
reference <- vapply(nom_vars, function(v) {
  if (betas[v] > 0) ranges[[v]][1] else ranges[[v]][2]
}, numeric(1))
max_logit_contribution <- vapply(nom_vars, function(v) {
  abs(betas[v]) * diff(ranges[[v]])
}, numeric(1))
logit_per_point <- max(max_logit_contribution) / 100
max_points <- max_logit_contribution / logit_per_point
total_max <- sum(max_points)
reference_logit <- unname(coef(final_model)["(Intercept)"] +
  sum(betas * reference))
nomogram <- data.frame(
  Predictor = nom_vars, Beta = unname(betas),
  Observed_min = vapply(ranges, `[`, numeric(1), 1),
  Observed_max = vapply(ranges, `[`, numeric(1), 2),
  Zero_point_reference = reference,
  Max_points = max_points,
  Logit_per_point = logit_per_point,
  Reference_logit = reference_logit)
write.csv(nomogram, file.path(out_dir, "06_nomogram_mapping_150.csv"),
          row.names = FALSE)
patient_points <- Reduce(`+`, lapply(nom_vars, function(v) {
  betas[v] * (primary[[v]] - reference[v]) / logit_per_point
}))
reconstructed <- plogis(reference_logit + patient_points * logit_per_point)
stopifnot(max(abs(reconstructed - apparent_probs)) < 1e-10)

nom_points <- function(v, value) {
  betas[v] * (value - reference[v]) / logit_per_point
}
nom_axis <- function(y, length_x, positions, labels, label,
                     point_cex = .83) {
  segments(0, y, length_x, y, col = "#31424D", lwd = 1.15)
  segments(positions, y, positions, y - .10,
           col = "#31424D", lwd = .95)
  for (i in seq_along(positions)) {
    align <- if (positions[i] < .025) 0 else
      if (positions[i] > .975) 1 else .5
    text(positions[i], y - .23, labels = labels[i],
         adj = c(align, .5), xpd = NA,
         cex = point_cex, col = "#26333B")
  }
  text(-.035, y + .01, labels = label, adj = 1,
       font = 2, cex = .88, xpd = NA)
}

draw_figure6 <- function() {
  # Two calibration/decision panels above one full-width nomogram. This keeps
  # the long variable labels inside panel C and legible at journal page size.
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE),
         widths = c(1, 1), heights = c(1.14, .86))
  par(bg = "white", family = "sans", las = 1)

  # A. Apparent observed rates (deciles) and the corrected logistic curve.
  par(mar = c(3.8, 4.0, 3.0, 1.0), mgp = c(2.3,.55,0),
      tcl = -.025, xaxs = "i", yaxs = "i")
  plot(NA_real_, xlim = c(0, 1), ylim = c(-.015, 1.02),
       xlab = "Predicted probability of ALNM",
       ylab = "Observed ALNM proportion",
       axes = FALSE, cex.lab = 1.08)
  axis(1, at = seq(0, 1, .2), cex.axis = .92)
  axis(2, at = seq(0, 1, .2), cex.axis = .92)
  box(col = "#88949B", lwd = .8)
  abline(0, 1, lty = 2, col = "#9EA8AE", lwd = 1.15)
  grid_p <- seq(.005, .995, length.out = 250)
  lines(grid_p, plogis(corrected_intercept +
                         corrected_slope * qlogis(grid_p)),
        col = "#21669B", lwd = 2.35)
  points(calibration_bins$Mean_predicted_risk,
         calibration_bins$Observed_ALNM_rate,
         pch = 21, bg = "white", col = "#D17C34", lwd = 1.25, cex = 1.12)
  legend("topleft", bty = "n", cex = .85,
         legend = c("Ideal", "Corrected logistic", "Observed deciles"),
         col = c("#9EA8AE", "#21669B", "#D17C34"),
         lty = c(2,1,NA), pch = c(NA,NA,21),
         pt.bg = c(NA,NA,"white"),
         lwd = c(1.15,2.35,1.25), y.intersp = 1.05)
  title(main = "A  Bootstrap validation", adj = 0,
        cex.main = 1.04, font.main = 2, line = 1.0)
  mtext(sprintf("n = 150; 67 events; 1,000 resamples"),
        side = 3, adj = 0, line = .08, cex = .75, col = "#5B656C")
  text(.96, .06, adj = c(1,0), cex = .79,
       labels = sprintf("AUC %.3f → %.3f", apparent_auc, corrected_auc),
       col = "#33424A")

  # B. Decision curve analysis. Extreme treat-all negative values are clipped
  # for display; numerical values in the CSV remain unchanged.
  par(mar = c(3.8, 3.8, 3.0, 1.0), mgp = c(2.3,.55,0),
      tcl = -.025, xaxs = "i", yaxs = "i")
  plot(NA_real_, xlim = c(0, .95), ylim = c(-.11, .50),
       xlab = "Threshold probability", ylab = "Net benefit",
       axes = FALSE, cex.lab = 1.08)
  axis(1, at = seq(0, .9, .1), cex.axis = .88)
  axis(2, at = seq(-.1, .5, .1), cex.axis = .9)
  box(col = "#88949B", lwd = .8)
  lines(thresholds, dca$Combined_model, col = "#1F5A82", lwd = 2.2)
  lines(thresholds, dca$Node_RADS_only, col = "#5285B3",
        lty = 2, lwd = 1.8)
  lines(thresholds, dca$Treat_all, col = "#999999",
        lty = 3, lwd = 1.35)
  abline(h = 0, col = "#333333", lwd = 1.1)
  legend("topright", bty = "n", cex = .85,
         legend = c("Combined model", "Node-RADS only",
                    "Treat all", "Treat none"),
         col = c("#1F5A82", "#5285B3", "#999999", "#333333"),
         lty = c(1,2,3,1), lwd = c(2.2,1.8,1.35,1.1),
         y.intersp = 1.06, seg.len = 1.9)
  title(main = "B  Decision curve analysis", adj = 0,
        cex.main = 1.04, font.main = 2, line = 1.0)
  mtext("Development cohort, n = 150", side = 3,
        adj = 0, line = .08, cex = .75, col = "#5B656C")

  # C. Four-predictor nomogram. Predictor rows use a common 0–100 point ruler;
  # total points and risk rows each use their own explicitly labeled scale.
  par(mar = c(1.2, 1.6, 3.0, .8), xaxs = "i", yaxs = "i")
  plot(NA_real_, xlim = c(-.25, 1.03), ylim = c(.52, 7.5),
       axes = FALSE, xlab = "", ylab = "", frame.plot = FALSE)
  nom_axis(6.9, 1, seq(0,1,.25), c("0","25","50","75","100"),
           "Points", point_cex = .86)
  for (v in nom_vars) {
    row <- switch(v, node_final = 6.0, size_mm = 5.02,
                  cADC = 4.04, peri_adc_new = 3.06)
    tick_values <- switch(v,
      node_final = 1:5,
      size_mm = c(7, 25, 45, 65, 85),
      cADC = seq(ranges[[v]][2], ranges[[v]][1], length.out = 5),
      peri_adc_new = c(ranges[[v]][1], mean(ranges[[v]]), ranges[[v]][2]))
    tick_pos <- nom_points(v, tick_values) / 100
    axis_end <- max_points[v] / 100
    labels <- switch(v,
      node_final = as.character(tick_values),
      size_mm = as.character(tick_values),
      cADC = sprintf("%.2f", tick_values),
      peri_adc_new = sprintf("%.2f", tick_values))
    label <- switch(v,
      node_final = "Node-RADS",
      size_mm = "Tumor diameter (mm)",
      cADC = "cADC",
      peri_adc_new = "Peritumoral ADC")
    nom_axis(row, axis_end, tick_pos, labels, label,
             point_cex = if (v == "node_final") .76 else .81)
  }
  total_ticks <- c(0,50,100,150,200,round(total_max))
  nom_axis(2.02, 1, total_ticks / total_max,
           as.character(total_ticks), "Total points", point_cex = .82)
  risk_ticks <- c(.05,.20,.50,.80,.95)
  risk_positions <- (qlogis(risk_ticks) - reference_logit) /
    logit_per_point / total_max
  keep <- risk_positions >= 0 & risk_positions <= 1
  nom_axis(1.06, 1, risk_positions[keep],
           sprintf("%.2f", risk_ticks[keep]),
           "ALNM probability", point_cex = .77)
  title(main = "C  Four-variable nomogram", adj = 0,
        cex.main = 1.04, font.main = 2, line = 1.0)
  mtext("Points reproduce the fitted logistic equation", side = 3,
        adj = 0, line = .08, cex = .75, col = "#5B656C")
}

png(file.path(out_dir, "Figure_6_validation_DCA_nomogram_150.png"),
    width = 12.5, height = 10.2, units = "in", res = 300,
    pointsize = 18, type = "cairo-png")
draw_figure6()
dev.off()

tiff(file.path(out_dir, "Figure_6_validation_DCA_nomogram_150.tiff"),
     width = 12.5, height = 10.2, units = "in", res = 300,
     pointsize = 18, compression = "lzw", type = "cairo")
draw_figure6()
dev.off()

pdf(file.path(out_dir, "Figure_6_validation_DCA_nomogram_150.pdf"),
    width = 12.5, height = 10.2, pointsize = 18, useDingbats = FALSE)
draw_figure6()
dev.off()

cat("\nFIGURE 6 RESULTS: 150 PATIENTS\n")
print(validation_summary, row.names = FALSE, digits = 7)
cat("DCA: combined model higher than all three comparators",
    "at every 0.01 grid point from 0.10 to 0.80.\n")
cat("Nomogram maximum total points:", total_max,
    "; reconstruction max error:",
    max(abs(reconstructed - apparent_probs)), "\n")
