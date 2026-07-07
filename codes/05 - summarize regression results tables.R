#!/usr/bin/env Rscript

# Create a table-only TeX summary from existing regression output CSV files.
# This script does not refit models.  It reads completed analysis outputs and
# writes G-regression_results_tables.tex.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

project_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)

result_root <- Sys.getenv(
  "PTB_SPLINE_RESULT_ROOT",
  unset = "C:/Users/moham/Box/Research_URES_ptb/analysis/06 - Splines"
)

out_tex <- Sys.getenv(
  "PTB_TABLE_TEX_OUT",
  unset = file.path(project_root, "G-regression_results_tables.tex")
)

fmt_num <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(
    is.na(x),
    "--",
    ifelse(abs(x) < 0.0005, "0.000", formatC(x, format = "f", digits = digits))
  )
}

fmt_p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "--", ifelse(x == 0, "0.000", formatC(x, format = "f", digits = 3)))
}

read_result <- function(...) {
  path <- file.path(result_root, ...)
  if (!file.exists(path)) stop("Missing result file: ", path)
  read_csv(path, show_col_types = FALSE)
}

table_lines <- function(rows) paste(rows, collapse = "\n")

delivery_table <- function() {
  global <- read_result("Test_results", "qr_globalTime_thetaPlusGamma_clusterBoot.csv")
  outcome <- read_result("Test_results", "qr_outcomeLevelPlusShape_clusterBoot.csv")
  shape <- read_result("Test_results", "qr_shapeOnly_clusterBoot.csv")

  rows <- global %>%
    left_join(outcome %>% select(var, tau, p_group = p_boot), by = c("var", "tau")) %>%
    left_join(shape %>% select(var, tau, p_shape = p_boot), by = c("var", "tau")) %>%
    transmute(
      line = paste0(
        var, " & ", fmt_num(tau, 2), " & ", fmt_p(p_boot), " & ",
        fmt_p(p_group), " & ", fmt_p(p_shape), " \\\\"
      )
    ) %>%
    pull(line)

  c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Delivery-aligned spline quantile regression tests. Time is \\(GA-GA_{\\mathrm{delivery}}\\). Entries are subject-level cluster-bootstrap \\(p\\)-values.}",
    "\\label{tab:delivery_spline_summary}",
    "\\scriptsize",
    "\\begin{tabular}{llccc}",
    "\\toprule",
    "Biomarker & \\(\\tau\\) & Global time & Level + shape & Shape only \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
}

ga_table <- function() {
  global <- read_result("Test_results_GA", "qr_constantOverTime_clusterBoot_GA.csv")
  group <- read_result("Test_results_GA", "qr_noGroupEffect_clusterBoot_GA.csv")

  rows <- global %>%
    left_join(group %>% select(var, tau, p_group = p_boot), by = c("var", "tau")) %>%
    transmute(
      line = paste0(var, " & ", fmt_num(tau, 2), " & ", fmt_p(p_boot), " & ", fmt_p(p_group), " \\\\")
    ) %>%
    pull(line)

  c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Gestational-age spline quantile regression tests using both groups. Entries are subject-level cluster-bootstrap \\(p\\)-values.}",
    "\\label{tab:ga_spline_summary}",
    "\\scriptsize",
    "\\begin{tabular}{llcc}",
    "\\toprule",
    "Biomarker & \\(\\tau\\) & Global time & No group effect \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
}

ftb_ga_table <- function() {
  ftb <- read_result("Test_results_GA_FTB", "qr_FTB_constantOverTime_clusterBoot_GA.csv")
  rows <- ftb %>%
    transmute(
      line = paste0(var, " & ", fmt_num(tau, 2), " & ", fmt_num(W_obs, 2), " & ", fmt_p(p_boot), " \\\\")
    ) %>%
    pull(line)

  c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{FTB-only gestational-age spline quantile regression tests for time variation. Entries are Wald statistics and subject-level cluster-bootstrap \\(p\\)-values.}",
    "\\label{tab:ftb_ga_spline_summary}",
    "\\scriptsize",
    "\\begin{tabular}{llcc}",
    "\\toprule",
    "Biomarker & \\(\\tau\\) & Wald statistic & \\(p_{\\mathrm{boot}}\\) \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
}

ald_phase_table <- function() {
  ald_dir <- file.path(result_root, "Test_results", "ALD_result")
  files <- list.files(ald_dir, pattern = "^ald_tshift_.*\\.csv$", full.names = TRUE)
  if (length(files) == 0) stop("No ALD phase-shift files found in ", ald_dir)

  dat <- bind_rows(lapply(files, read_csv, show_col_types = FALSE)) %>%
    mutate(tau_num = as.numeric(tau)) %>%
    arrange(factor(biomarker, levels = c("SWS", "Midband", "Intercept", "AC", "Slope")), tau_num)

  rows <- dat %>%
    transmute(
      line = paste0(
        biomarker, " & ", fmt_num(tau, 2), " & ", fmt_num(phi_hat), " & (",
        fmt_num(phi_ci_low), ", ", fmt_num(phi_ci_high), ") & ",
        fmt_p(p_shift), " & ", fmt_num(sig2_hat), " & ", fmt_p(p_het), " \\\\"
      )
    ) %>%
    pull(line)

  c(
    "\\begin{table*}[!htbp]",
    "\\centering",
    "\\caption{PTB-versus-FTB ALD--Student-\\(t\\) random-shift phase model. All rows use 1082 observations from 521 subjects.}",
    "\\label{tab:ald_phase_shift_summary}",
    "\\scriptsize",
    "\\begin{tabular}{llccccc}",
    "\\toprule",
    "Biomarker & \\(\\tau\\) & \\(\\widehat\\phi_\\tau\\) & 95\\% CI for \\(\\phi_\\tau\\) & \\(p_{\\mathrm{shift}}\\) & \\(\\widehat\\sigma_{\\phi,\\tau}^2\\) & \\(p_{\\mathrm{het}}\\) \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table*}"
  )
}

ftb_ald_table <- function() {
  dat <- read_result("Test_results", "ALD_FTB_only", "summary_FTB_random_shift_only.csv") %>%
    mutate(tau_num = as.numeric(tau)) %>%
    arrange(factor(biomarker, levels = c("SWS", "Midband", "Intercept", "AC", "Slope")), tau_num)

  rows <- dat %>%
    transmute(
      line = paste0(
        biomarker, " & ", fmt_num(tau, 2), " & ", fmt_num(sig2_hat), " & (",
        fmt_num(sig2_ci_low), ", ", fmt_num(sig2_ci_high), ") & ", fmt_p(p_het), " \\\\"
      )
    ) %>%
    pull(line)

  c(
    "\\begin{table*}[!htbp]",
    "\\centering",
    "\\caption{FTB-only ALD--Student-\\(t\\) random-shift variance model. All rows use 939 observations from 454 FTB subjects.}",
    "\\label{tab:ftb_ald_phase_variance_summary}",
    "\\scriptsize",
    "\\begin{tabular}{llccc}",
    "\\toprule",
    "Biomarker & \\(\\tau\\) & \\(\\widehat\\sigma_{\\phi,\\tau}^2\\) & 95\\% CI for \\(\\sigma_{\\phi,\\tau}^2\\) & \\(p_{\\mathrm{het}}\\) \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table*}"
  )
}

semiparametric_table <- function() {
  path <- Sys.getenv(
    "PTB_SEMIPARAM_CSV",
    unset = file.path(project_root, "Test_results_semiparametric_FTB", "ftb_semiparametric_phase_variance_all.csv")
  )
  if (!file.exists(path)) stop("Missing semiparametric result file: ", path)

  dat <- read_csv(path, show_col_types = FALSE) %>%
    mutate(
      tau_num = as.numeric(tau),
      biomarker = factor(biomarker, levels = c("SWS", "Midband", "Intercept", "AC", "Slope")),
      method = factor(method, levels = c("quadratic", "spline")),
      ci = paste0("(", fmt_num(ci_low), ", ", fmt_num(ci_high), ")")
    ) %>%
    arrange(method, tau_num, biomarker)

  n_obs <- unique(dat$n)
  n_id <- unique(dat$n_id)
  b_ok <- unique(dat$B_ok)
  if (length(n_obs) != 1 || length(n_id) != 1) {
    n_text <- "available FTB observations"
  } else {
    n_text <- paste0(n_obs, " observations from ", n_id, " FTB subjects")
  }
  if (length(b_ok) == 1) {
    boot_text <- paste0(" Each row uses ", b_ok, " subject-level bootstrap replicates.")
  } else {
    boot_text <- " Bootstrap replicate counts are reported by row."
  }

  rows <- dat %>%
    transmute(
      line = paste0(
        as.character(biomarker), " & ", fmt_num(tau, 2), " & ", as.character(method),
        " & ", fmt_num(sigma_b2), " & ", fmt_num(sigma_b2_positive),
        " & ", ci, " & ", fmt_p(p_two_sided_normal), " & ",
        fmt_p(p_upper_boot), " & ", B_ok, " \\\\"
      )
    ) %>%
    pull(line)

  c(
    "\\begin{table*}[!htbp]",
    "\\centering",
    paste0(
      "\\caption{FTB-only semiparametric phase-variance estimates for all biomarkers. All rows use ",
      n_text,
      ".",
      boot_text,
      "}"
    ),
    "\\label{tab:ftb_semiparametric_phase_variance_summary}",
    "\\scriptsize",
    "\\begin{tabular}{lllcccccc}",
    "\\toprule",
    "Biomarker & \\(\\tau\\) & Method & \\(\\widehat\\sigma_b^2\\) & \\(\\widehat\\sigma_{b,+}^2\\) & 95\\% CI & Normal \\(p\\) & \\(p_{\\le 0}^{\\mathrm{boot}}\\) & \\(B_{ok}\\) \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table*}"
  )
}

interpretation_key <- function() {
  c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Abbreviated interpretation key.}",
    "\\label{tab:regression_interpretation_key}",
    "\\scriptsize",
    "\\begin{tabular}{lll}",
    "\\toprule",
    "Block & Main result pattern & Tabular conclusion \\\\",
    "\\midrule",
    "Delivery-aligned splines & Global time significant, shape mostly not significant & Time structure without stable shape separation \\\\",
    "Gestational-age splines & Global time significant, group mostly not significant & Gestational trajectories vary over time, weak group separation \\\\",
    "PTB-versus-FTB phase shift & Average shift mainly SWS; heterogeneity often positive & Timing signal is biomarker-specific and model-based \\\\",
    "FTB-only ALD variance & Positive model-based random-shift variance & FTB timing heterogeneity under ALD model \\\\",
    "FTB-only semiparametric & Quadratic CIs include zero; spline raw estimates mostly nonpositive & No robust distribution-light evidence for positive phase variance \\\\",
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
}

all_lines <- c(
  "\\section{G. Regression Results Tables}",
  "",
  delivery_table(),
  "",
  ga_table(),
  "",
  ftb_ga_table(),
  "",
  ald_phase_table(),
  "",
  ftb_ald_table(),
  "",
  semiparametric_table(),
  "",
  interpretation_key()
)

writeLines(all_lines, out_tex)
message("Wrote ", out_tex)
