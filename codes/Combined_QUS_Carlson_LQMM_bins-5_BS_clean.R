## Combined QUS + Carlson LQMM analysis
##
## Standalone Rscript version of Combined_QUS_Carlson_LQMM_bins-5_BS.Rmd.
## The script reads every git/data/V1_V7_longformat_k*.csv file by default
## and writes all outputs under:
##   git/results/Combined_QUS_Carlson_LQMM_bins-5_BS

required_packages <- c("readr", "dplyr", "ggplot2", "splines", "lqmm")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Install required R package(s) before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(splines)
  library(lqmm)
})

analysis_name <- "Combined_QUS_Carlson_LQMM_bins-5_BS"

illini_orange <- "#FF5F05"
illini_blue <- "#13294B"
illini_white <- "#FFFFFF"
storm <- "#707372"
harvest <- "#FCB316"
prairie <- "#006230"
source_alpha <- 0.72

get_script_dir <- function() {
  file_arg <- "--file="
  cmd_args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub(file_arg, "", cmd_args[startsWith(cmd_args, file_arg)])

  if (length(script_path) > 0 && nzchar(script_path[[1]])) {
    return(dirname(normalizePath(script_path[[1]], mustWork = TRUE)))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(active_path)) {
      return(dirname(normalizePath(active_path, mustWork = TRUE)))
    }
  }

  normalizePath(getwd(), mustWork = TRUE)
}

check_required_columns <- function(data, required, data_name) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(
      data_name,
      " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

parse_k_values <- function(k_text) {
  if (!nzchar(k_text)) return(integer(0))
  as.integer(trimws(unlist(strsplit(k_text, ","))))
}

discover_qus_files <- function(data_dir, k_values = integer(0)) {
  qus_files <- list.files(
    data_dir,
    pattern = "^V1_V7_longformat_k[0-9]+\\.csv$",
    full.names = TRUE
  )

  if (length(qus_files) == 0) {
    stop("No QUS files found matching V1_V7_longformat_k*.csv in ", data_dir, call. = FALSE)
  }

  qus_config <- data.frame(
    k = as.integer(sub("^V1_V7_longformat_k([0-9]+)\\.csv$", "\\1", basename(qus_files))),
    path = qus_files,
    stringsAsFactors = FALSE
  ) %>%
    arrange(k)

  if (length(k_values) > 0) {
    qus_config <- qus_config %>% filter(k %in% k_values)
  }

  if (nrow(qus_config) == 0) {
    stop("No QUS files matched requested PTB_K_VALUES: ", paste(k_values, collapse = ", "), call. = FALSE)
  }

  qus_config
}

prepare_qus <- function(qus_path) {
  qus_raw <- read_csv(qus_path, show_col_types = FALSE)
  check_required_columns(qus_raw, c("Participant_ID", "GA", "SWS", "mPTB", "Outcome"), "QUS data")

  qus_raw %>%
    filter(as.character(mPTB) != "1") %>%
    filter(as.character(Outcome) == "0") %>%
    transmute(
      Participant_ID = paste0("QUS_", Participant_ID),
      GA = as.numeric(GA) * 39,
      Avg_SWS = as.numeric(SWS),
      source = "QUS"
    ) %>%
    filter(is.finite(GA), is.finite(Avg_SWS), !is.na(Participant_ID))
}

prepare_carlson <- function(carlson_path, collapse_per_exam = TRUE) {
  carlson_raw <- read.table(carlson_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  check_required_columns(carlson_raw, c("ID", "Exam", "GA", "P", "MP1", "MP2", "M"), "Carlson data")

  carlson <- carlson_raw %>%
    mutate(Avg_SWS = (P + MP1 + MP2 + M) / 4)

  if (isTRUE(collapse_per_exam)) {
    carlson <- carlson %>%
      group_by(ID, Exam, GA) %>%
      summarise(Avg_SWS = mean(Avg_SWS, na.rm = TRUE), .groups = "drop")
  }

  carlson %>%
    transmute(
      Participant_ID = paste0("CARL_", ID),
      GA = as.numeric(GA),
      Avg_SWS = as.numeric(Avg_SWS),
      source = "Carlson"
    ) %>%
    filter(is.finite(GA), is.finite(Avg_SWS), !is.na(Participant_ID))
}

fit_lqmm_curves <- function(data, response, basis_names, grid_ga, ga_ns, tau_values) {
  fixed_form <- as.formula(paste(response, "~", paste(basis_names, collapse = " + ")))
  grid_x <- cbind(Intercept = 1, predict(ga_ns, grid_ga))

  fits <- vector("list", length(tau_values))
  names(fits) <- paste0("tau_", tau_values)

  predictions <- bind_rows(lapply(seq_along(tau_values), function(i) {
    tau <- tau_values[[i]]

    fit <- lqmm(
      fixed = fixed_form,
      random = ~ 1,
      group = Participant_ID,
      tau = tau,
      data = data,
      control = lqmmControl(LP_max_iter = 1000)
    )

    fits[[i]] <<- fit
    beta <- as.numeric(coef(fit))

    data.frame(
      GA = grid_ga,
      fit = as.numeric(grid_x %*% beta),
      tau = factor(paste0("tau= ", tau), levels = paste0("tau= ", tau_values)),
      stringsAsFactors = FALSE
    )
  }))

  coefficients <- bind_rows(lapply(seq_along(fits), function(i) {
    beta <- coef(fits[[i]])
    data.frame(
      tau = tau_values[[i]],
      term = names(beta),
      estimate = as.numeric(beta),
      stringsAsFactors = FALSE
    )
  }))

  list(fits = fits, predictions = predictions, coefficients = coefficients)
}

save_png_plot <- function(plot, basename, plot_dir, width = 9, height = 6) {
  ggsave(
    file.path(plot_dir, paste0(basename, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}

make_linear_plot <- function(combined, predictions, bin_plot_positions, bin_labels, y_top) {
  ggplot(combined, aes(x = GA_mid, y = Avg_SWS, group = GA_mid)) +
    geom_boxplot(outlier.alpha = 0.12, fill = storm, color = storm, alpha = source_alpha, width = 3) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 95,
      size = 11,
      color = illini_white,
      show.legend = FALSE
    ) +
    geom_line(
      data = predictions,
      aes(x = GA, y = fit, color = tau, linetype = tau, group = tau),
      inherit.aes = FALSE,
      linewidth = 1.3
    ) +
    scale_color_manual(values = tau_colors) +
    scale_linetype_manual(values = tau_linetypes) +
    scale_x_continuous(breaks = bin_plot_positions, labels = bin_labels) +
    labs(
      title = "Avg_SWS vs GA (Combined QUS + Carlson)",
      subtitle = "5-week boxplots with LQMM population quantiles",
      x = "Gestational age (weeks)",
      y = "Avg_SWS",
      color = "Population quantiles",
      linetype = "Population quantiles"
    ) +
    coord_cartesian(xlim = c(5, 100), ylim = c(0, y_top)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

make_log_plot <- function(combined_log, predictions, bin_plot_positions, bin_labels, y_low, y_high) {
  ggplot(combined_log, aes(x = GA_mid, y = log_SWS, group = GA_mid)) +
    geom_boxplot(outlier.alpha = 0.12, fill = storm, color = storm, alpha = source_alpha, width = 3) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 95,
      size = 11,
      color = illini_white,
      show.legend = FALSE
    ) +
    geom_line(
      data = predictions,
      aes(x = GA, y = fit, color = tau, linetype = tau, group = tau),
      inherit.aes = FALSE,
      linewidth = 1.3
    ) +
    scale_color_manual(values = tau_colors) +
    scale_linetype_manual(values = tau_linetypes) +
    scale_x_continuous(breaks = bin_plot_positions, labels = bin_labels) +
    labs(
      title = "log(Avg_SWS) vs GA (Combined QUS + Carlson)",
      subtitle = "5-week boxplots with LQMM population quantiles",
      x = "Gestational age (weeks)",
      y = "log(Avg_SWS)",
      color = "Population quantiles",
      linetype = "Population quantiles"
    ) +
    coord_cartesian(xlim = c(5, 100), ylim = c(y_low, y_high)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

make_source_plot <- function(combined, predictions, bin_plot_positions, bin_labels, y_top) {
  ggplot(
    combined,
    aes(x = GA_mid, y = Avg_SWS, group = interaction(GA_mid, source), fill = source)
  ) +
    geom_boxplot(
      outlier.alpha = 0.12,
      color = storm,
      alpha = source_alpha,
      width = 3,
      position = position_dodge(width = 2)
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 95,
      size = 8,
      color = illini_white,
      position = position_dodge(width = 2),
      show.legend = FALSE
    ) +
    geom_line(
      data = predictions,
      aes(x = GA, y = fit, color = tau, linetype = tau, group = tau),
      inherit.aes = FALSE,
      linewidth = 1.3
    ) +
    scale_fill_manual(values = source_colors) +
    scale_color_manual(values = tau_colors) +
    scale_linetype_manual(values = tau_linetypes) +
    scale_x_continuous(breaks = bin_plot_positions, labels = bin_labels) +
    labs(
      title = "Avg_SWS vs GA by source (combined LQMM curves overlaid)",
      x = "Gestational age (weeks)",
      y = "Avg_SWS",
      fill = "Source",
      color = "Population quantiles",
      linetype = "Population quantiles"
    ) +
    coord_cartesian(xlim = c(6, 96), ylim = c(0, y_top)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

make_source_log_plot <- function(combined_log, predictions, bin_plot_positions, bin_labels, y_low, y_high) {
  ggplot(
    combined_log,
    aes(x = GA_mid, y = log_SWS, group = interaction(GA_mid, source), fill = source)
  ) +
    geom_boxplot(
      outlier.alpha = 0.12,
      color = storm,
      alpha = source_alpha,
      width = 3,
      position = position_dodge(width = 2)
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 95,
      size = 8,
      color = illini_white,
      position = position_dodge(width = 2),
      show.legend = FALSE
    ) +
    geom_line(
      data = predictions,
      aes(x = GA, y = fit, color = tau, linetype = tau, group = tau),
      inherit.aes = FALSE,
      linewidth = 1.3
    ) +
    scale_fill_manual(values = source_colors) +
    scale_color_manual(values = tau_colors) +
    scale_linetype_manual(values = tau_linetypes) +
    scale_x_continuous(breaks = bin_plot_positions, labels = bin_labels) +
    labs(
      title = "log(Avg_SWS) vs GA by source (combined LQMM curves overlaid)",
      x = "Gestational age (weeks)",
      y = "log(Avg_SWS)",
      fill = "Source",
      color = "Population quantiles",
      linetype = "Population quantiles"
    ) +
    coord_cartesian(xlim = c(8, 100), ylim = c(y_low, y_high)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

run_analysis_for_k <- function(k, qus_path, carlson, config) {
  suffix <- paste0("_k", k)
  message("Running k=", k, " using QUS data: ", qus_path)

  qus <- prepare_qus(qus_path)

  combined <- bind_rows(qus, carlson) %>%
    mutate(Participant_ID = as.factor(Participant_ID)) %>%
    arrange(Participant_ID, GA)

  ga_ns <- ns(combined$GA, knots = config$spline_knots)
  basis_df <- as.data.frame(ga_ns)
  names(basis_df) <- paste0("ns", seq_len(ncol(basis_df)))
  basis_names <- names(basis_df)

  combined <- bind_cols(combined, basis_df) %>%
    mutate(
      GA_bin = cut(GA, breaks = config$bin_breaks, right = FALSE, include.lowest = TRUE),
      GA_mid = config$bin_breaks[as.integer(GA_bin)] + config$bin_position_offset
    ) %>%
    filter(!is.na(GA_mid))

  data_summary <- combined %>%
    group_by(source) %>%
    summarise(
      k = k,
      n_rows = n(),
      n_participants = n_distinct(Participant_ID),
      GA_min = min(GA, na.rm = TRUE),
      GA_max = max(GA, na.rm = TRUE),
      Avg_SWS_min = min(Avg_SWS, na.rm = TRUE),
      Avg_SWS_median = median(Avg_SWS, na.rm = TRUE),
      Avg_SWS_max = max(Avg_SWS, na.rm = TRUE),
      .groups = "drop"
    )

  bin_counts <- combined %>%
    group_by(source, GA_bin, GA_mid) %>%
    summarise(
      k = k,
      n_rows = n(),
      n_participants = n_distinct(Participant_ID),
      .groups = "drop"
    ) %>%
    arrange(GA_mid, source)

  write_csv(data_summary, file.path(config$table_dir, paste0("data_summary", suffix, ".csv")))
  write_csv(bin_counts, file.path(config$table_dir, paste0("bin_counts", suffix, ".csv")))
  write_csv(combined, file.path(config$table_dir, paste0("combined_analysis_data", suffix, ".csv")))

  grid_ga <- seq(min(combined$GA), max(combined$GA_mid), length.out = 300)
  linear_results <- fit_lqmm_curves(
    data = combined,
    response = "Avg_SWS",
    basis_names = basis_names,
    grid_ga = grid_ga,
    ga_ns = ga_ns,
    tau_values = config$tau_values
  )

  combined_log <- combined %>%
    filter(Avg_SWS > 0) %>%
    mutate(log_SWS = log(Avg_SWS))

  grid_ga_log <- seq(min(combined_log$GA), max(combined_log$GA_mid), length.out = 300)
  log_results <- fit_lqmm_curves(
    data = combined_log,
    response = "log_SWS",
    basis_names = basis_names,
    grid_ga = grid_ga_log,
    ga_ns = ga_ns,
    tau_values = config$tau_values
  )

  write_csv(linear_results$predictions, file.path(config$table_dir, paste0("lqmm_predictions_linear", suffix, ".csv")))
  write_csv(log_results$predictions, file.path(config$table_dir, paste0("lqmm_predictions_log", suffix, ".csv")))
  write_csv(linear_results$coefficients, file.path(config$table_dir, paste0("lqmm_coefficients_linear", suffix, ".csv")))
  write_csv(log_results$coefficients, file.path(config$table_dir, paste0("lqmm_coefficients_log", suffix, ".csv")))

  saveRDS(linear_results$fits, file.path(config$model_dir, paste0("lqmm_models_linear", suffix, ".rds")))
  saveRDS(log_results$fits, file.path(config$model_dir, paste0("lqmm_models_log", suffix, ".rds")))

  y_top <- quantile(combined$Avg_SWS, 0.995, na.rm = TRUE) * 1.05
  y_log_low <- min(combined_log$log_SWS, na.rm = TRUE)
  y_log_high <- quantile(combined_log$log_SWS, 0.995, na.rm = TRUE)

  combined_plot <- make_linear_plot(
    combined,
    linear_results$predictions,
    config$bin_plot_positions,
    config$bin_labels,
    y_top
  )
  combined_plot_log <- make_log_plot(
    combined_log,
    log_results$predictions,
    config$bin_plot_positions,
    config$bin_labels,
    y_log_low,
    y_log_high
  )
  source_plot <- make_source_plot(
    combined,
    linear_results$predictions,
    config$bin_plot_positions,
    config$bin_labels,
    y_top
  )
  source_plot_log <- make_source_log_plot(
    combined_log,
    log_results$predictions,
    config$bin_plot_positions,
    config$bin_labels,
    y_log_low,
    y_log_high
  )

  save_png_plot(combined_plot, paste0("combined_lqmm_linear", suffix), config$plot_dir)
  save_png_plot(combined_plot_log, paste0("combined_lqmm_log", suffix), config$plot_dir)
  save_png_plot(source_plot, paste0("combined_lqmm_linear_by_source", suffix), config$plot_dir)
  save_png_plot(source_plot_log, paste0("combined_lqmm_log_by_source", suffix), config$plot_dir)

  manifest <- data.frame(
    field = c(
      "analysis_name",
      "run_time",
      "k",
      "project_root",
      "qus_path",
      "carlson_path",
      "output_dir",
      "collapse_carlson_per_exam",
      "tau_values",
      "spline_knots",
      "bin_width",
      "bin_position_offset",
      "illini_orange",
      "illini_blue",
      "carlson_data_color",
      "qus_data_color",
      "source_alpha",
      "quantile_tau_0_2_color",
      "quantile_tau_0_5_color",
      "quantile_tau_0_8_color",
      "box_median_color",
      "plot_format"
    ),
    value = c(
      analysis_name,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      as.character(k),
      config$project_root,
      normalizePath(qus_path, mustWork = TRUE),
      config$carlson_path,
      normalizePath(config$output_dir, mustWork = FALSE),
      as.character(config$collapse_carlson_per_exam),
      paste(config$tau_values, collapse = ", "),
      paste(config$spline_knots, collapse = ", "),
      as.character(config$bin_width),
      as.character(config$bin_position_offset),
      illini_orange,
      illini_blue,
      source_colors[["Carlson"]],
      source_colors[["QUS"]],
      as.character(source_alpha),
      tau_colors[["tau= 0.2"]],
      tau_colors[["tau= 0.5"]],
      tau_colors[["tau= 0.8"]],
      illini_white,
      "png"
    ),
    stringsAsFactors = FALSE
  )
  write_csv(manifest, file.path(config$table_dir, paste0("run_manifest", suffix, ".csv")))

  data_summary
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
data_dir <- normalizePath(Sys.getenv("PTB_DATA_DIR", unset = file.path(project_root, "data")), mustWork = TRUE)
results_root <- Sys.getenv("PTB_RESULTS_DIR", unset = file.path(project_root, "results"))
output_dir <- Sys.getenv(
  "PTB_COMBINED_QUS_CARLSON_OUT_DIR",
  unset = file.path(results_root, analysis_name)
)

plot_dir <- file.path(output_dir, "plots")
table_dir <- file.path(output_dir, "tables")
model_dir <- file.path(output_dir, "models")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

carlson_path <- Sys.getenv("PTB_CARLSON_DATA_PATH", unset = file.path(data_dir, "carlson_et_al_data.txt"))
if (!file.exists(carlson_path)) stop("Missing Carlson data file: ", carlson_path, call. = FALSE)

collapse_carlson_per_exam <- TRUE
tau_values <- c(0.2, 0.5, 0.8)
tau_colors <- c(
  "tau= 0.2" = prairie,
  "tau= 0.5" = illini_orange,
  "tau= 0.8" = illini_blue
)
tau_linetypes <- c(
  "tau= 0.2" = "solid",
  "tau= 0.5" = "solid",
  "tau= 0.8" = "solid"
)
source_colors <- c("Carlson" = storm, "QUS" = harvest)
spline_knots <- c(16, 20, 24, 31, 36, 40, 44, 48, 64)
bin_width <- 5
bin_position_offset <- 2
bin_breaks <- seq(6, 96, by = bin_width)
bin_plot_positions <- head(bin_breaks, -1) + bin_position_offset
bin_labels <- paste0("[", head(bin_breaks, -1), ",", tail(bin_breaks, -1), ")")
bin_labels[length(bin_labels)] <- sub("\\)$", "]", bin_labels[length(bin_labels)])

k_values <- parse_k_values(Sys.getenv("PTB_K_VALUES", unset = ""))
qus_config <- discover_qus_files(data_dir, k_values)

message("Reading Carlson data: ", carlson_path)
carlson <- prepare_carlson(carlson_path, collapse_carlson_per_exam)

config <- list(
  project_root = project_root,
  output_dir = output_dir,
  plot_dir = plot_dir,
  table_dir = table_dir,
  model_dir = model_dir,
  carlson_path = normalizePath(carlson_path, mustWork = TRUE),
  collapse_carlson_per_exam = collapse_carlson_per_exam,
  tau_values = tau_values,
  spline_knots = spline_knots,
  bin_width = bin_width,
  bin_position_offset = bin_position_offset,
  bin_breaks = bin_breaks,
  bin_plot_positions = bin_plot_positions,
  bin_labels = bin_labels
)

message("QUS k values: ", paste(qus_config$k, collapse = ", "))
summaries <- bind_rows(lapply(seq_len(nrow(qus_config)), function(i) {
  run_analysis_for_k(
    k = qus_config$k[[i]],
    qus_path = qus_config$path[[i]],
    carlson = carlson,
    config = config
  )
}))

write_csv(summaries, file.path(table_dir, "data_summary_all_k.csv"))

message("Done. PNG outputs written to: ", normalizePath(plot_dir, mustWork = FALSE))
message("Tables and models written to: ", normalizePath(output_dir, mustWork = FALSE))
