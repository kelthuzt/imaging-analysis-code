# Node-RADS / ADC breast MRI study: reproducible analysis
# Part 01: import and primary cohort audit (additional parts will be added in order).
# This script contains no patient names or individual data.
#
# For submission: put node数据更新--9-8.xlsx beside this script, install readxl,
# and run: Rscript NodeRADS_analysis.R
# If the workbook is stored elsewhere, set NODE_INPUT_XLSX to its path.
# The optional NODE_INPUT_CSV environment variable is used only for an exact,
# deidentified cell export when the local R runtime lacks readxl.

options(stringsAsFactors = FALSE)
source_xlsx <- Sys.getenv("NODE_INPUT_XLSX", unset = "node数据更新--9-8.xlsx")
source_csv <- Sys.getenv("NODE_INPUT_CSV", unset = "")
output_dir <- Sys.getenv("NODE_OUTPUT_DIR", unset = "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

to_number <- function(x) suppressWarnings(as.numeric(trimws(as.character(x))))
to_label <- function(x) {
  y <- trimws(as.character(x))
  y[is.na(x) | !nzchar(y)] <- NA_character_
  y
}

read_study_cells <- function(xlsx, csv = "") {
  if (nzchar(csv)) {
    if (!file.exists(csv)) stop("CSV fallback not found: ", csv)
    dat <- read.csv(csv, check.names = FALSE, na.strings = c("", "NA"),
                    fileEncoding = "UTF-8")
  } else {
    if (!file.exists(xlsx)) stop("Workbook not found: ", xlsx)
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Install the readxl package to import the original XLSX.")
    }
    raw <- as.data.frame(readxl::read_excel(xlsx, sheet = 1,
                                            .name_repair = "minimal"))
    # Source workbook columns C:AB, selected by their original Excel position.
    indices <- c(3, 4, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28)
    dat <- raw[, indices, drop = FALSE]
    names(dat) <- c("lvi", "menopausal", "bmi", "age", "laterality",
                    "path_type_code", "path_type", "er", "pr", "her2",
                    "ki67", "subtype", "alnm", "node_senior",
                    "node_junior", "node_final", "tic", "tumor_adc",
                    "ipsi_adc", "contra_adc", "size_raw", "peri_adc_new",
                    "peri_adc_old", "size_mm")
    dat$source_row <- seq_len(nrow(dat)) + 1L
  }
  expected <- c("source_row", "lvi", "laterality", "path_type", "her2",
                "subtype", "alnm", "node_senior", "node_junior",
                "node_final", "tumor_adc", "ipsi_adc", "contra_adc",
                "peri_adc_new", "size_mm")
  if (!all(expected %in% names(dat))) {
    stop("Input is missing columns: ",
         paste(setdiff(expected, names(dat)), collapse = ", "))
  }
  numeric_cols <- c("source_row", "lvi", "menopausal", "bmi", "age",
                    "laterality", "path_type_code", "er", "pr", "her2",
                    "ki67", "alnm", "node_senior", "node_junior",
                    "node_final", "tic", "tumor_adc", "ipsi_adc",
                    "contra_adc", "peri_adc_new", "peri_adc_old", "size_mm")
  for (column in numeric_cols) dat[[column]] <- to_number(dat[[column]])
  dat$subtype <- to_label(dat$subtype)
  dat$path_type <- to_label(dat$path_type)
  dat
}

all_records <- read_study_cells(source_xlsx, source_csv)

# The original workbook is already a partly screened extract. It cannot
# independently verify the 195 screened patients or earlier exclusion counts.
# The following rules act on every patient record present in the supplied XLSX.
all_records$flag_bilateral <- is.na(all_records$laterality) |
  !(all_records$laterality %in% c(1, 2))
all_records$flag_in_situ <- is.na(all_records$path_type) |
  grepl("原位癌|in situ", all_records$path_type, ignore.case = TRUE)
all_records$flag_missing_outcome <- is.na(all_records$alnm) |
  !(all_records$alnm %in% c(0, 1))
all_records$flag_unmeasurable_contralateral_adc <-
  is.na(all_records$contra_adc) | !is.finite(all_records$contra_adc) |
  all_records$contra_adc <= 0
other_required <- c("node_senior", "node_junior", "node_final", "tumor_adc",
                    "ipsi_adc", "peri_adc_new", "size_mm")
all_records$flag_other_imaging_missing <- apply(
  all_records[, other_required, drop = FALSE], 1L,
  function(x) any(is.na(x) | !is.finite(x) | x <= 0)
)
all_records$included_primary <- !(
  all_records$flag_bilateral | all_records$flag_in_situ |
  all_records$flag_missing_outcome |
  all_records$flag_unmeasurable_contralateral_adc |
  all_records$flag_other_imaging_missing
)

primary <- all_records[all_records$included_primary, , drop = FALSE]
primary$cADC <- primary$ipsi_adc / primary$contra_adc
primary$rADC <- primary$ipsi_adc / primary$tumor_adc
stopifnot(all(is.finite(primary$cADC)), all(is.finite(primary$rADC)))

# Molecular subtype, HER2, and postoperative LVI are NOT primary-cohort filters.
# Raw ipsilateral, contralateral, and tumor ADC are retained for descriptive
# analyses and ratio calculation, but excluded from the multivariable candidate
# predictor set at the author's request. Use peri_adc_new, not peri_adc_old.
model_adc_variables <- c("cADC", "rADC", "peri_adc_new")
excluded_raw_adc_predictors <- c("ipsi_adc", "contra_adc", "tumor_adc")

audit <- all_records[, c("source_row", "alnm", "included_primary",
                         "flag_bilateral", "flag_in_situ",
                         "flag_missing_outcome",
                         "flag_unmeasurable_contralateral_adc",
                         "flag_other_imaging_missing"), drop = FALSE]
audit$subtype_recorded <- !is.na(all_records$subtype)
audit$lvi_recorded <- !is.na(all_records$lvi)
write.csv(audit, file.path(output_dir, "01_cohort_row_audit.csv"),
          row.names = FALSE, na = "")

screened_from_manuscript <- 195L
counts <- data.frame(
  metric = c("Screened in manuscript (requires source verification)",
             "Records in supplied workbook", "Bilateral within workbook",
             "In situ or missing tumor histology within workbook",
             "Missing coded axillary outcome within workbook",
             "Unmeasurable contralateral ADC within workbook",
             "Other required imaging missing within workbook",
             "PRIMARY COHORT", "ALNM positive", "ALNM negative",
             "Female (author-confirmed from source records)",
             "Subtype label available in primary cohort",
             "LVI available in primary cohort"),
  n = c(screened_from_manuscript, nrow(all_records),
        sum(all_records$flag_bilateral), sum(all_records$flag_in_situ),
        sum(all_records$flag_missing_outcome),
        sum(all_records$flag_unmeasurable_contralateral_adc),
        sum(all_records$flag_other_imaging_missing), nrow(primary),
        sum(primary$alnm == 1), sum(primary$alnm == 0),
        150L,
        sum(!is.na(primary$subtype)), sum(!is.na(primary$lvi)))
)
if (nrow(primary) != 150L) {
  stop("Recheck author-confirmed sex count if the main cohort is no longer 150.")
}
write.csv(counts, file.path(output_dir, "01_cohort_counts.csv"),
          row.names = FALSE)
# Earlier exclusions were not present as rows in the supplied workbook.
# Their numbers below come from the manuscript, with the author's confirmed
# correction from 15 to 12 missing postoperative axillary pathology results.
flow <- data.frame(
  reason = c("Ductal carcinoma in situ",
             "Missing postoperative axillary pathology",
             "Axillary ADC coverage absent",
             "Contralateral reference-node ADC unmeasurable",
             "Bilateral breast cancer", "MRI performed after surgery",
             "Nondiagnostic MRI quality"),
  n = c(17L, 12L, 4L, 8L, 2L, 1L, 1L)
)
stopifnot(sum(flow$n) == screened_from_manuscript - nrow(primary))
write.csv(flow, file.path(output_dir, "01_flow_counts_author_reconciled.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "01_R_sessionInfo.txt"))
print(counts, row.names = FALSE)
cat("Source Excel rows excluded from the primary cohort:",
    paste(all_records$source_row[!all_records$included_primary],
          collapse = ", "), "\n")
cat("Source Excel rows retained despite missing subtype:",
    paste(primary$source_row[is.na(primary$subtype)], collapse = ", "), "\n")

# Part 02: Table 1 and the opening paragraph of Results section 3.1.
# Descriptive/test rules match the manuscript: assess normality separately
# within the ALNM-negative and -positive groups. Use Welch's independent t
# test and mean (SD) if BOTH groups have Shapiro-Wilk P >= 0.05; otherwise
# use the Wilcoxon rank-sum test and median (IQR). Categorical variables use
# Pearson chi-square without continuity correction if all expected counts
# are >= 5 and Fisher's exact test otherwise. The Node-RADS score is ordinal.
format_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) "<0.001" else sprintf("%.3f", p)
}
mean_sd <- function(x) sprintf("%.2f ± %.2f", mean(x), sd(x))
median_iqr <- function(x) {
  q <- quantile(x, c(.25, .5, .75), na.rm = TRUE, names = FALSE)
  sprintf("%.2f (%.2f, %.2f)", q[2], q[1], q[3])
}
table1_rows <- list()
add_table1_row <- function(label, display_all, display0, display1, p, test,
                           n_all, n0, n1,
                           shapiro0 = NA_real_, shapiro1 = NA_real_) {
  table1_rows[[length(table1_rows) + 1L]] <<- data.frame(
    Characteristic = label, Overall = display_all,
    `Non-ALNM` = display0, ALNM = display1,
    `N overall` = n_all, `N non-ALNM` = n0, `N ALNM` = n1,
    P = format_p(p), Test = test,
    `Shapiro P non-ALNM` = shapiro0, `Shapiro P ALNM` = shapiro1,
    check.names = FALSE
  )
}
add_continuous <- function(variable, label, ordinal = FALSE) {
  x0 <- primary[[variable]][primary$alnm == 0]
  x1 <- primary[[variable]][primary$alnm == 1]
  x0 <- x0[is.finite(x0)]; x1 <- x1[is.finite(x1)]
  p0 <- if (length(x0) >= 3L) shapiro.test(x0)$p.value else NA_real_
  p1 <- if (length(x1) >= 3L) shapiro.test(x1)$p.value else NA_real_
  is_normal <- !ordinal && is.finite(p0) && is.finite(p1) &&
    p0 >= 0.05 && p1 >= 0.05
  if (is_normal) {
    test <- "Welch independent-samples t test"
    p <- t.test(x0, x1, var.equal = FALSE)$p.value
    display_all <- mean_sd(c(x0, x1))
    display0 <- mean_sd(x0); display1 <- mean_sd(x1)
  } else {
    test <- "Wilcoxon rank-sum test"
    p <- wilcox.test(x0, x1, exact = FALSE)$p.value
    display_all <- median_iqr(c(x0, x1))
    display0 <- median_iqr(x0); display1 <- median_iqr(x1)
  }
  add_table1_row(label, display_all, display0, display1, p, test,
                 length(x0) + length(x1), length(x0), length(x1), p0, p1)
}
add_categorical <- function(variable, label, levels, labels) {
  values <- factor(primary[[variable]], levels = levels, labels = labels)
  outcome <- factor(primary$alnm, levels = c(0, 1))
  tab <- table(values, outcome)
  expected <- suppressWarnings(chisq.test(tab, correct = FALSE)$expected)
  if (any(expected < 5)) {
    test <- "Fisher exact test"
    p <- fisher.test(tab)$p.value
  } else {
    test <- "Pearson chi-square test"
    p <- chisq.test(tab, correct = FALSE)$p.value
  }
  totals <- colSums(tab)
  add_table1_row(label, "", "", "", p, test,
                 sum(totals), totals[1], totals[2])
  for (i in seq_along(labels)) {
    value <- function(n, total) sprintf("%d (%.1f%%)", n, 100 * n / total)
    add_table1_row(paste0("  ", labels[i]),
                   value(sum(tab[i, ]), sum(totals)),
                   value(tab[i, 1], totals[1]), value(tab[i, 2], totals[2]),
                   NA_real_, "", sum(totals), totals[1], totals[2])
  }
}

add_continuous("age", "Age, years")
# Sex is absent from the supplied analysis workbook. The author checked the
# clinical source records and confirmed that all 150 included patients are
# female. It is reported descriptively; no group comparison is possible.
add_table1_row("Female sex, n (%)", "150 (100.0%)", "83 (100.0%)",
               "67 (100.0%)", NA_real_, "Not applicable (no variation)",
               150L, 83L, 67L)
add_continuous("bmi", "BMI, kg/m²")
add_categorical("menopausal", "Menopausal status", c(0, 1),
                c("Premenopausal", "Postmenopausal"))
add_categorical("laterality", "Side", c(1, 2), c("Left", "Right"))
add_categorical("path_type_code", "Histological type", c(1, 2),
                c("Non-special type", "Special type"))
add_categorical("er", "ER status", c(1, 2), c("Positive", "Negative"))
add_categorical("pr", "PR status", c(1, 2), c("Positive", "Negative"))
add_categorical("her2", "HER2 status", c(1, 2), c("Positive", "Negative"))
add_categorical("ki67", "Ki-67", c(1, 2), c("<20%", "≥20%"))
add_categorical("tic", "TIC type", c(1, 2, 3),
                c("Washout", "Plateau", "Persistent"))
add_continuous("size_mm", "Maximum tumor diameter, mm")
add_continuous("tumor_adc", "Tumor ADC, 10^-3 mm²/s")
add_continuous("ipsi_adc", "Ipsilateral node ADC, 10^-3 mm²/s")
add_continuous("contra_adc", "Contralateral node ADC, 10^-3 mm²/s")
add_continuous("cADC", "cADC")
add_continuous("rADC", "rADC")
add_continuous("peri_adc_new", "Peritumoral ADC (new), 10^-3 mm²/s")
add_continuous("node_final", "Final Node-RADS score", ordinal = TRUE)

table1 <- do.call(rbind, table1_rows)
write.csv(table1, file.path(output_dir, "02_Table1_all_variables.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")
cat("\nTABLE 1, 150 patients (all original main-table variables):\n")
print(table1[, c("Characteristic", "Overall", "Non-ALNM", "ALNM", "P")],
      row.names = FALSE, right = FALSE)
cat("\nData checks for Table 1:\n")
cat("All 150 patients were female, as verified by the author against",
    "source records; sex was not a column in the supplied spreadsheet.\n")
cat("HER2 available:", sum(!is.na(primary$her2)), "of", nrow(primary),
    "; among ALNM-negative:",
    sum(!is.na(primary$her2) & primary$alnm == 0),
    "; among ALNM-positive:",
    sum(!is.na(primary$her2) & primary$alnm == 1), "\n")

# Part 03: Results section 3.2, patient-level Node-RADS reader assessment.
# The same 150 patients are used throughout. The original reader scores are
# kept separate from the adjudicated/final scores. The ordinal ROC uses the
# empirical Mann-Whitney AUC and DeLong's analytic 95% CI (including ties).
# Binary paired comparisons use exact McNemar P values, calculated as exact
# binomial tests of the discordant pairs. No optional ROC/kappa package is
# required; the formulae below use base R and have been checked against the
# published 147-patient final-score AUC/CI before updating to 150 patients.
stopifnot(nrow(primary) == 150L,
          sum(primary$alnm == 1) == 67L,
          sum(primary$alnm == 0) == 83L)
score_variables <- c("node_final", "node_senior", "node_junior")
for (variable in score_variables) {
  x <- primary[[variable]]
  stopifnot(!anyNA(x), all(x %in% 1:5))
}

auc_delong <- function(outcome, score) {
  positive <- score[outcome == 1]
  negative <- score[outcome == 0]
  comparison <- outer(positive, negative, `>`) +
    0.5 * outer(positive, negative, `==`)
  case_contribution <- rowMeans(comparison)
  control_contribution <- colMeans(comparison)
  estimate <- mean(comparison)
  standard_error <- sqrt(var(case_contribution) / length(positive) +
                         var(control_contribution) / length(negative))
  c(auc = estimate,
    ci_lower = max(0, estimate - qnorm(0.975) * standard_error),
    ci_upper = min(1, estimate + qnorm(0.975) * standard_error),
    se = standard_error)
}

reader_performance <- function(outcome, score, reader) {
  thresholds <- 1:5
  sensitivity <- vapply(thresholds,
                        function(t) mean(score[outcome == 1] >= t), numeric(1))
  specificity <- vapply(thresholds,
                        function(t) mean(score[outcome == 0] < t), numeric(1))
  cutoff <- thresholds[which.max(sensitivity + specificity - 1)]
  predicted <- score >= cutoff
  tp <- sum(predicted & outcome == 1)
  fn <- sum(!predicted & outcome == 1)
  tn <- sum(!predicted & outcome == 0)
  fp <- sum(predicted & outcome == 0)
  roc <- auc_delong(outcome, score)
  data.frame(Reader = reader, N = length(outcome),
             ALNM_positive = sum(outcome == 1),
             ALNM_negative = sum(outcome == 0),
             AUC = unname(roc["auc"]),
             CI_lower = unname(roc["ci_lower"]),
             CI_upper = unname(roc["ci_upper"]),
             Youden_cutoff_ge = cutoff,
             TP = tp, FN = fn, TN = tn, FP = fp,
             Accuracy = (tp + tn) / length(outcome),
             Sensitivity = tp / (tp + fn),
             Specificity = tn / (tn + fp))
}

reader_results <- do.call(rbind, list(
  reader_performance(primary$alnm, primary$node_final, "Final/adjudicated"),
  reader_performance(primary$alnm, primary$node_senior, "Senior original"),
  reader_performance(primary$alnm, primary$node_junior, "Junior original")
))
stopifnot(all(reader_results$Youden_cutoff_ge == 3L))
write.csv(reader_results, file.path(output_dir, "03_NodeRADS_performance.csv"),
          row.names = FALSE)

# Quadratic-weighted Cohen kappa: disagreement weights (i-j)^2 / (5-1)^2.
score_levels <- 1:5
pair_counts <- table(Senior = factor(primary$node_senior, levels = score_levels),
                     Junior = factor(primary$node_junior, levels = score_levels))
observed <- pair_counts / sum(pair_counts)
expected <- outer(rowSums(observed), colSums(observed))
disagreement_weight <- outer(score_levels, score_levels,
                             function(i, j) ((i - j) / 4)^2)
weighted_kappa <- 1 - sum(disagreement_weight * observed) /
  sum(disagreement_weight * expected)
pair_output <- as.data.frame(pair_counts, responseName = "N")
write.csv(pair_output, file.path(output_dir, "03_NodeRADS_score_pairs.csv"),
          row.names = FALSE)

distribution <- do.call(rbind, lapply(score_variables, function(variable) {
  counts <- as.data.frame(table(
    Score = factor(primary[[variable]], levels = score_levels),
    ALNM = factor(primary$alnm, levels = c(0, 1))), responseName = "N")
  data.frame(Reader = variable, counts)
}))
write.csv(distribution, file.path(output_dir,
                                  "03_NodeRADS_score_distribution.csv"),
          row.names = FALSE)

exact_mcnemar <- function(senior_binary, junior_binary) {
  stopifnot(length(senior_binary) == length(junior_binary),
            !anyNA(senior_binary), !anyNA(junior_binary))
  senior_only <- sum(senior_binary & !junior_binary)
  junior_only <- sum(!senior_binary & junior_binary)
  p <- if (senior_only + junior_only == 0L) 1 else
    binom.test(senior_only, senior_only + junior_only,
               p = 0.5, alternative = "two.sided")$p.value
  c(senior_only = senior_only, junior_only = junior_only, P = p)
}

senior_positive <- primary$node_senior >= 3
junior_positive <- primary$node_junior >= 3
senior_correct <- senior_positive == (primary$alnm == 1)
junior_correct <- junior_positive == (primary$alnm == 1)
positive_cases <- primary$alnm == 1
negative_cases <- primary$alnm == 0
score_discordance <- primary$node_senior != primary$node_junior
threshold_discordance <- senior_positive != junior_positive

accuracy_pair <- exact_mcnemar(senior_correct, junior_correct)
detection_pair <- exact_mcnemar(senior_positive[positive_cases],
                                junior_positive[positive_cases])
score_one_pair <- exact_mcnemar(primary$node_senior[negative_cases] == 1,
                                primary$node_junior[negative_cases] == 1)
threshold_pair <- exact_mcnemar(senior_correct[threshold_discordance],
                                junior_correct[threshold_discordance])
discordance_by_outcome <- table(
  Discrepant = factor(score_discordance, levels = c(FALSE, TRUE)),
  ALNM = factor(primary$alnm, levels = c(0, 1)))
stopifnot(all(chisq.test(discordance_by_outcome,
                        correct = TRUE)$expected >= 5))
# Yates correction matches the manuscript's original 147-patient P = 0.002.
discordance_p <- chisq.test(discordance_by_outcome,
                            correct = TRUE)$p.value

comparisons <- data.frame(
  Metric = c("Quadratic-weighted Cohen kappa",
             "Accuracy: senior vs junior, exact McNemar",
             "Detection among ALNM-positive, exact McNemar",
             "Score 1 among ALNM-negative, exact McNemar",
             "All-score disagreements",
             "Cross-threshold disagreements",
             "Correct among cross-threshold, exact McNemar",
             "ALNM in discrepant vs concordant, Yates-corrected chi-square",
             "Jointly missed ALNM",
             "Jointly detected ALNM",
             "Final score correct among cross-threshold cases"),
  Value = c(weighted_kappa,
            sum(senior_correct) / nrow(primary),
            sum(senior_positive[positive_cases]) / sum(positive_cases),
            sum(primary$node_senior[negative_cases] == 1) /
              sum(negative_cases),
            sum(score_discordance), sum(threshold_discordance),
            sum(senior_correct[threshold_discordance]),
            sum(primary$alnm[score_discordance] == 1) /
              sum(score_discordance),
            sum(positive_cases & !senior_positive & !junior_positive),
            sum(positive_cases & senior_positive & junior_positive),
            sum((primary$node_final[threshold_discordance] >= 3) ==
                  (primary$alnm[threshold_discordance] == 1))),
  Comparator = c(NA_real_,
                 sum(junior_correct) / nrow(primary),
                 sum(junior_positive[positive_cases]) / sum(positive_cases),
                 sum(primary$node_junior[negative_cases] == 1) /
                   sum(negative_cases),
                 NA_real_, NA_real_,
                 sum(junior_correct[threshold_discordance]),
                 sum(primary$alnm[!score_discordance] == 1) /
                   sum(!score_discordance),
                 NA_real_, NA_real_, NA_real_),
  N = c(nrow(primary), nrow(primary), sum(positive_cases),
        sum(negative_cases), nrow(primary), nrow(primary),
        sum(threshold_discordance), nrow(primary),
        sum(positive_cases), sum(positive_cases),
        sum(threshold_discordance)),
  Senior_only = c(NA_real_, accuracy_pair["senior_only"],
                  detection_pair["senior_only"],
                  score_one_pair["senior_only"], NA_real_, NA_real_,
                  threshold_pair["senior_only"], NA_real_, NA_real_,
                  NA_real_, NA_real_),
  Junior_only = c(NA_real_, accuracy_pair["junior_only"],
                  detection_pair["junior_only"],
                  score_one_pair["junior_only"], NA_real_, NA_real_,
                  threshold_pair["junior_only"], NA_real_, NA_real_,
                  NA_real_, NA_real_),
  P = c(NA_real_, accuracy_pair["P"], detection_pair["P"],
        score_one_pair["P"], NA_real_, NA_real_, threshold_pair["P"],
        discordance_p, NA_real_, NA_real_, NA_real_))
write.csv(comparisons, file.path(output_dir, "03_NodeRADS_comparisons.csv"),
          row.names = FALSE, na = "")

# Prespecified seven imaging measures in the positive cases jointly missed
# versus jointly detected. Raw ipsilateral/contralateral/tumor ADC can be used
# for descriptive comparisons; the author excludes these raw values only from
# multivariable model candidates. Holm correction covers all seven measures.
missed <- positive_cases & !senior_positive & !junior_positive
detected <- positive_cases & senior_positive & junior_positive
imaging_measures <- c("ipsi_adc", "contra_adc", "tumor_adc", "cADC",
                      "rADC", "peri_adc_new", "size_mm")
missed_results <- do.call(rbind, lapply(imaging_measures, function(variable) {
  a <- primary[[variable]][missed]
  b <- primary[[variable]][detected]
  normal <- shapiro.test(a)$p.value >= 0.05 &&
    shapiro.test(b)$p.value >= 0.05
  if (normal) {
    p <- t.test(a, b, var.equal = FALSE)$p.value
    test <- "Welch independent-samples t test"
  } else {
    p <- wilcox.test(a, b, exact = FALSE)$p.value
    test <- "Wilcoxon rank-sum test"
  }
  data.frame(Measure = variable, Missed_N = length(a),
             Detected_N = length(b),
             Missed_median = median(a), Detected_median = median(b),
             Missed_mean = mean(a), Detected_mean = mean(b),
             Test = test, P = p)
}))
missed_results$Holm_P <- p.adjust(missed_results$P, method = "holm")
write.csv(missed_results,
          file.path(output_dir, "03_NodeRADS_jointly_missed.csv"),
          row.names = FALSE)

cat("\nSECTION 3.2, UPDATED 150-PATIENT NODE-RADS RESULTS:\n")
print(reader_results, row.names = FALSE, digits = 5)
print(comparisons, row.names = FALSE, digits = 5)
print(missed_results, row.names = FALSE, digits = 5)

# Part 04: Results 3.3 and Table 2: preoperative imaging model (n = 150).
# Candidates are prespecified imaging variables, including the current/new
# peritumoral ADC column. Raw tumor, ipsilateral-node, and contralateral-node
# ADC are descriptive only and cannot re-enter this candidate set. TIC is
# screened as a categorical imaging variable. Pathology/subtype/LVI are not
# primary-model candidates.
stopifnot(nrow(primary) == 150L,
          !anyNA(primary[, c("alnm", "node_final", "size_mm", "cADC",
                               "rADC", "peri_adc_new", "tic")]))
candidate_imaging <- c("cADC", "rADC", "peri_adc_new", "size_mm",
                       "node_final", "tic")
stopifnot(!any(candidate_imaging %in% excluded_raw_adc_predictors))

# For continuous variables, use two-sided Wald P and Wald 95% CIs; TIC is
# evaluated as a whole using the likelihood-ratio test against an intercept.
# An OR is shown for each 0.1 increment in the three ADC measures, for each
# 1-mm increment in tumor diameter, and for each 1-point Node-RADS increment.
or_scale <- c(cADC = 0.1, rADC = 0.1, peri_adc_new = 0.1,
              size_mm = 1, node_final = 1)
variable_label <- c(cADC = "cADC (per 0.1)",
                    rADC = "rADC (per 0.1)",
                    peri_adc_new = "Peritumoral ADC (per 0.1 × 10⁻³ mm²/s)",
                    size_mm = "Maximum tumor diameter (per 1 mm)",
                    node_final = "Node-RADS (per 1 point)")
univ_fits <- setNames(lapply(candidate_imaging, function(v) {
  if (v == "tic") {
    glm(alnm ~ factor(tic), data = primary, family = binomial())
  } else {
    glm(reformulate(v, response = "alnm"), data = primary,
        family = binomial())
  }
}), candidate_imaging)
univ_p <- vapply(candidate_imaging, function(v) {
  fit <- univ_fits[[v]]
  if (v == "tic") {
    anova(glm(alnm ~ 1, data = primary, family = binomial()),
          fit, test = "Chisq")$`Pr(>Chi)`[2L]
  } else {
    coef(summary(fit))[v, "Pr(>|z|)"]
  }
}, numeric(1))
screened <- setdiff(candidate_imaging[univ_p < 0.10], "tic")
# If categorical TIC passes P < .10 in a future dataset, stop rather than
# quietly omit it: factor terms require a block-level backward selection.
if ("tic" %in% candidate_imaging[univ_p < 0.10]) {
  stop("TIC passed screening: extend backward selection to handle its factor block")
}

current_terms <- screened
selection_trace <- data.frame(Stage = character(), Removed = character(),
                              Largest_Wald_P = numeric())
repeat {
  fit <- glm(reformulate(current_terms, response = "alnm"),
             data = primary, family = binomial())
  term_p <- coef(summary(fit))[current_terms, "Pr(>|z|)"]
  if (all(term_p < 0.05)) break
  remove_term <- names(which.max(term_p))
  selection_trace <- rbind(selection_trace,
                           data.frame(Stage = paste(current_terms,
                                                    collapse = "+"),
                                      Removed = remove_term,
                                      Largest_Wald_P = max(term_p)))
  current_terms <- setdiff(current_terms, remove_term)
  if (length(current_terms) == 0L) stop("No imaging predictor remained")
}
final_model <- fit
stopifnot(setequal(current_terms,
                   c("cADC", "peri_adc_new", "size_mm", "node_final")))
baseline_model <- glm(alnm ~ node_final + cADC, data = primary,
                      family = binomial())

wald_or <- function(fit, variable, increment) {
  row <- coef(summary(fit))[variable, ]
  beta <- row[["Estimate"]] * increment
  se <- row[["Std. Error"]] * increment
  c(OR = exp(beta), CI_lower = exp(beta - qnorm(0.975) * se),
    CI_upper = exp(beta + qnorm(0.975) * se),
    P = row[["Pr(>|z|)"]])
}
table2 <- do.call(rbind, lapply(names(or_scale), function(variable) {
  univariate <- wald_or(univ_fits[[variable]], variable,
                        or_scale[[variable]])
  multivariate <- if (variable %in% current_terms) {
    wald_or(final_model, variable, or_scale[[variable]])
  } else rep(NA_real_, 4L)
  display_or <- function(x) {
    if (anyNA(x)) "—" else sprintf("%.3f (%.3f–%.3f)",
                                      x["OR"], x["CI_lower"], x["CI_upper"])
  }
  data.frame(Factor = variable_label[[variable]],
             Univariable_OR_95CI = display_or(univariate),
             Univariable_P = format_p(univariate["P"]),
             Multivariable_OR_95CI = display_or(multivariate),
             Multivariable_P = if (is.na(multivariate["P"])) "—" else
               format_p(multivariate["P"]),
             Univariate_OR = univariate["OR"],
             Univariate_CI_lower = univariate["CI_lower"],
             Univariate_CI_upper = univariate["CI_upper"],
             Univariate_P_exact = univariate["P"],
             Multivariate_OR = multivariate["OR"],
             Multivariate_CI_lower = multivariate["CI_lower"],
             Multivariate_CI_upper = multivariate["CI_upper"],
             Multivariate_P_exact = multivariate["P"],
             check.names = FALSE)
}))
write.csv(table2, file.path(output_dir, "04_Table2_logistic_150.csv"),
          row.names = FALSE, na = "")
write.csv(data.frame(Candidate = candidate_imaging,
                     Univariable_P = univ_p,
                     Passed_P_below_0_10 = univ_p < 0.10),
          file.path(output_dir, "04_model_candidate_screen_150.csv"),
          row.names = FALSE)
write.csv(selection_trace,
          file.path(output_dir, "04_model_backward_selection_150.csv"),
          row.names = FALSE)

# Two-sided paired DeLong comparison for ROC curves on the same patients;
# case/control placement vectors incorporate ties. The AUC difference
# percentile CI uses 1,000 paired stratified case/control bootstrap samples
# of the original fitted scores (no model refitting). This CI quantifies
# sampling variation in apparent discrimination, not internal validation.
paired_auc_delong <- function(outcome, first, second) {
  pos <- which(outcome == 1); neg <- which(outcome == 0)
  placement <- function(score) {
    cmp <- outer(score[pos], score[neg], `>`) +
      0.5 * outer(score[pos], score[neg], `==`)
    list(auc = mean(cmp), v10 = rowMeans(cmp), v01 = colMeans(cmp))
  }
  x <- placement(first); y <- placement(second)
  cov_auc <- cov(x$v10, y$v10) / length(pos) +
    cov(x$v01, y$v01) / length(neg)
  var_diff <- var(x$v10) / length(pos) + var(x$v01) / length(neg) +
    var(y$v10) / length(pos) + var(y$v01) / length(neg) - 2 * cov_auc
  difference <- x$auc - y$auc
  c(first_auc = x$auc, second_auc = y$auc,
    difference = difference, z = difference / sqrt(var_diff),
    P = 2 * pnorm(-abs(difference / sqrt(var_diff))))
}
auc_bootstrap_difference <- function(outcome, first, second, B = 1000L,
                                     seed = 20260928L) {
  set.seed(seed)
  pos <- which(outcome == 1); neg <- which(outcome == 0)
  diffs <- replicate(B, {
    chosen <- c(sample(pos, length(pos), replace = TRUE),
                sample(neg, length(neg), replace = TRUE))
    auc_delong(outcome[chosen], first[chosen])["auc"] -
      auc_delong(outcome[chosen], second[chosen])["auc"]
  })
  unname(quantile(diffs, c(0.025, 0.975)))
}
youden_probability <- function(outcome, probabilities) {
  u <- sort(unique(probabilities))
  cutoffs <- c(-Inf, (head(u, -1) + tail(u, -1)) / 2, Inf)
  sensitivities <- vapply(cutoffs,
                          function(x) mean(probabilities[outcome == 1] >= x),
                          numeric(1))
  specificities <- vapply(cutoffs,
                          function(x) mean(probabilities[outcome == 0] < x),
                          numeric(1))
  best <- which.max(sensitivities + specificities - 1)
  pred <- probabilities >= cutoffs[best]
  c(Threshold = cutoffs[best], Sensitivity = sensitivities[best],
    Specificity = specificities[best],
    Accuracy = mean(pred == (outcome == 1)),
    TP = sum(pred & outcome == 1), FN = sum(!pred & outcome == 1),
    TN = sum(!pred & outcome == 0), FP = sum(pred & outcome == 0))
}

pred_final <- predict(final_model, type = "response")
pred_base <- predict(baseline_model, type = "response")
full_auc <- auc_delong(primary$alnm, pred_final)
base_auc <- auc_delong(primary$alnm, pred_base)
auc_comparison <- paired_auc_delong(primary$alnm, pred_final, pred_base)
auc_boot_ci <- auc_bootstrap_difference(primary$alnm, pred_final,
                                        pred_base)
lrt <- anova(baseline_model, final_model, test = "Chisq")
cutoff <- youden_probability(primary$alnm, pred_final)
performance <- data.frame(
  Model = c("Node-RADS + cADC baseline", "Four-variable combined"),
  N = 150L,
  AUC = c(base_auc["auc"], full_auc["auc"]),
  CI_lower = c(base_auc["ci_lower"], full_auc["ci_lower"]),
  CI_upper = c(base_auc["ci_upper"], full_auc["ci_upper"]),
  Cutoff_ge = c(NA_real_, cutoff["Threshold"]),
  Sensitivity = c(NA_real_, cutoff["Sensitivity"]),
  Specificity = c(NA_real_, cutoff["Specificity"]),
  Accuracy = c(NA_real_, cutoff["Accuracy"]),
  TP = c(NA_integer_, cutoff["TP"]),
  FN = c(NA_integer_, cutoff["FN"]),
  TN = c(NA_integer_, cutoff["TN"]),
  FP = c(NA_integer_, cutoff["FP"]), row.names = NULL)
comparison <- data.frame(
  Metric = c("AUC difference full minus baseline",
             "Paired two-sided DeLong P", "Bootstrap AUC difference CI lower",
             "Bootstrap AUC difference CI upper", "Likelihood ratio chi-square",
             "Likelihood ratio df", "Likelihood ratio P"),
  Value = c(auc_comparison["difference"], auc_comparison["P"],
            auc_boot_ci[1], auc_boot_ci[2],
            lrt$Deviance[2], lrt$Df[2], lrt$`Pr(>Chi)`[2]))
equation <- data.frame(Term = names(coef(final_model)),
                       Coefficient = unname(coef(final_model)))
write.csv(performance,
          file.path(output_dir, "04_model_performance_150.csv"),
          row.names = FALSE, na = "")
write.csv(comparison,
          file.path(output_dir, "04_model_comparisons_150.csv"),
          row.names = FALSE)
write.csv(equation,
          file.path(output_dir, "04_model_equation_150.csv"),
          row.names = FALSE)
cat("\nSECTION 3.3 MODEL, 150 PATIENTS:\n")
print(selection_trace, row.names = FALSE)
print(table2[, 1:5], row.names = FALSE, right = FALSE)
print(equation, row.names = FALSE, digits = 8)
print(performance, row.names = FALSE, digits = 8)
print(comparison, row.names = FALSE, digits = 8)
