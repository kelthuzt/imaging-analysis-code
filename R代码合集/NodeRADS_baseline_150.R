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
