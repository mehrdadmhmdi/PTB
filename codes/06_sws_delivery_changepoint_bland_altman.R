#!/usr/bin/env Rscript

# SWS analyses on the delivery-normalized clock:
#   1) near-delivery minimum SWS range and descriptive change point;
#   2) sPTB trajectories after normalizing to time-to-delivery;
#   3) Bland-Altman agreement of SWS across k-averaging values.

rm(list = ls())

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

ILLINI_ORANGE <- "#FF5F05"
ILLINI_BLUE <- "#13294B"

CLOSE_WINDOW_WEEKS <- as.numeric(Sys.getenv("PTB_CLOSE_WINDOW_WEEKS", "10"))
CHANGE_WINDOW_WEEKS <- as.numeric(Sys.getenv("PTB_CHANGE_WINDOW_WEEKS", "24"))
TIME_BIN_WIDTH_WEEKS <- as.numeric(Sys.getenv("PTB_TIME_BIN_WIDTH_WEEKS", "2"))
CHANGE_GRID_STEP_WEEKS <- as.numeric(Sys.getenv("PTB_CHANGE_GRID_STEP_WEEKS", "0.25"))
CHANGE_MIN_N <- as.integer(Sys.getenv("PTB_CHANGE_MIN_N", "30"))

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]),
                                 winslash = "/", mustWork = FALSE)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    if (!is.null(active_path) && nzchar(active_path)) {
      return(dirname(normalizePath(active_path, winslash = "/", mustWork = FALSE)))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

parse_k_values <- function(default) {
  value <- Sys.getenv("PTB_SWS_K_VALUES", unset = "")
  if (!nzchar(value)) return(default)

  if (grepl("^\\s*\\d+\\s*:\\s*\\d+\\s*$", value)) {
    bounds <- as.integer(strsplit(gsub("\\s+", "", value), ":", fixed = TRUE)[[1]])
    return(seq(bounds[[1]], bounds[[2]]))
  }

  values <- as.integer(strsplit(gsub("\\s+", "", value), ",", fixed = TRUE)[[1]])
  values <- values[!is.na(values)]
  if (length(values) == 0) stop("PTB_SWS_K_VALUES did not contain any valid integers.")
  sort(unique(values))
}

available_k_values <- function(data_dir) {
  files <- list.files(data_dir, pattern = "^V1_V7_longformat_k[0-9]+\\.csv$")
  values <- as.integer(sub("^V1_V7_longformat_k([0-9]+)\\.csv$", "\\1", files))
  sort(values[!is.na(values)])
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

ga_multiplier <- function(x) {
  mx <- suppressWarnings(max(as_num(x), na.rm = TRUE))
  if (is.finite(mx) && mx <= 2.5) 39 else 1
}

outcome_label <- function(x) {
  raw <- trimws(as.character(x))
  dplyr::case_when(
    raw %in% c("1", "PTB", "sPTB", "ptb", "sptb", "Preterm", "preterm") ~ "sPTB",
    raw %in% c("0", "FTB", "ftb", "Fullterm", "Full term", "Term", "term") ~ "FTB",
    TRUE ~ raw
  )
}

period_label <- function(x) {
  dplyr::case_when(
    x < 0 ~ "Prepartum",
    x > 0 ~ "Postpartum",
    TRUE ~ "Delivery"
  )
}

min_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else min(x)
}

max_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

q_or_na <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, probs = p, names = FALSE, na.rm = TRUE))
}

fmt_range <- function(x, digits = 2) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_character_)
  paste0(formatC(min(x), format = "f", digits = digits),
         " to ",
         formatC(max(x), format = "f", digits = digits))
}

add_time_bin <- function(dat, width) {
  dat %>%
    mutate(
      time_bin_left = floor(time_to_delivery_weeks / width) * width,
      time_bin_right = time_bin_left + width,
      time_bin_mid = 0.5 * (time_bin_left + time_bin_right),
      time_bin = sprintf("[%s,%s)",
                         formatC(time_bin_left, format = "f", digits = 1),
                         formatC(time_bin_right, format = "f", digits = 1))
    )
}

read_long_k <- function(k, data_dir) {
  path <- file.path(data_dir, sprintf("V1_V7_longformat_k%d.csv", k))
  if (!file.exists(path)) stop("Missing long-format file for k=", k, ": ", path)

  raw <- read_csv(path, show_col_types = FALSE)
  ga_mult <- ga_multiplier(raw$GA)
  delivery_mult <- ga_multiplier(raw$GA_delivery)

  raw %>%
    mutate(
      k = as.integer(k),
      Participant_ID = as.character(Participant_ID),
      Visit_ID = as.integer(Visit_ID),
      Outcome_label = outcome_label(Outcome),
      GA_weeks = as_num(GA) * ga_mult,
      GA_delivery_weeks = as_num(GA_delivery) * delivery_mult,
      time_to_delivery_weeks = GA_weeks - GA_delivery_weeks,
      period = factor(period_label(time_to_delivery_weeks),
                      levels = c("Prepartum", "Delivery", "Postpartum")),
      SWS = as_num(SWS)
    ) %>%
    filter(
      !is.na(Participant_ID),
      !is.na(Visit_ID),
      !is.na(Outcome_label),
      is.finite(GA_weeks),
      is.finite(GA_delivery_weeks),
      is.finite(time_to_delivery_weeks),
      is.finite(SWS)
    )
}

fit_piecewise_change <- function(dat, close_window, min_n, step) {
  d <- dat %>%
    transmute(x = time_to_delivery_weeks, y = SWS) %>%
    filter(is.finite(x), is.finite(y)) %>%
    arrange(x)

  n <- nrow(d)
  if (n < min_n || length(unique(round(d$x, 8))) < 6) {
    return(tibble(
      n_obs = n,
      n_unique_time = length(unique(round(d$x, 8))),
      change_point_weeks = NA_real_,
      nearest_observed_to_change_weeks = NA_real_,
      change_point_in_delivery_gap = NA,
      rss_linear = NA_real_,
      rss_piecewise = NA_real_,
      delta_aic_vs_linear = NA_real_,
      slope_before_change = NA_real_,
      slope_after_change = NA_real_,
      model_min_time_within_close_window = NA_real_,
      model_min_sws_within_close_window = NA_real_
    ))
  }

  x <- d$x
  y <- d$y
  lo <- max(min(x) + step, as.numeric(quantile(x, 0.10, names = FALSE)))
  hi <- min(max(x) - step, as.numeric(quantile(x, 0.90, names = FALSE)))
  candidates <- seq(lo, hi, by = step)

  if (length(candidates) == 0) {
    candidates <- sort(unique(x))
    candidates <- candidates[candidates > min(x) & candidates < max(x)]
  }
  if (length(candidates) == 0) {
    return(tibble(
      n_obs = n,
      n_unique_time = length(unique(round(d$x, 8))),
      change_point_weeks = NA_real_,
      nearest_observed_to_change_weeks = NA_real_,
      change_point_in_delivery_gap = NA,
      rss_linear = NA_real_,
      rss_piecewise = NA_real_,
      delta_aic_vs_linear = NA_real_,
      slope_before_change = NA_real_,
      slope_after_change = NA_real_,
      model_min_time_within_close_window = NA_real_,
      model_min_sws_within_close_window = NA_real_
    ))
  }

  linear_fit <- lm(y ~ x)
  rss_linear <- sum(residuals(linear_fit)^2)
  aic_linear <- n * log(max(rss_linear / n, .Machine$double.eps)) + 2 * 2

  best <- NULL
  for (cp in candidates) {
    h <- pmax(0, x - cp)
    fit <- tryCatch(lm(y ~ x + h), error = function(e) NULL)
    if (is.null(fit)) next
    rss <- sum(residuals(fit)^2)
    aic <- n * log(max(rss / n, .Machine$double.eps)) + 2 * length(coef(fit))
    if (is.null(best) || aic < best$aic) {
      best <- list(cp = cp, fit = fit, rss = rss, aic = aic)
    }
  }

  if (is.null(best)) {
    return(tibble(
      n_obs = n,
      n_unique_time = length(unique(round(d$x, 8))),
      change_point_weeks = NA_real_,
      nearest_observed_to_change_weeks = NA_real_,
      change_point_in_delivery_gap = NA,
      rss_linear = rss_linear,
      rss_piecewise = NA_real_,
      delta_aic_vs_linear = NA_real_,
      slope_before_change = NA_real_,
      slope_after_change = NA_real_,
      model_min_time_within_close_window = NA_real_,
      model_min_sws_within_close_window = NA_real_
    ))
  }

  beta <- coef(best$fit)
  slope_before <- unname(beta[["x"]])
  slope_after <- unname(beta[["x"]] + beta[["h"]])
  nearest_obs <- min(abs(x - best$cp))
  closest_pre <- max_or_na(x[x < 0])
  closest_post <- min_or_na(x[x > 0])
  in_gap <- is.finite(closest_pre) &&
    is.finite(closest_post) &&
    best$cp > closest_pre &&
    best$cp < closest_post

  grid_lo <- max(min(x), -close_window)
  grid_hi <- min(max(x), close_window)
  if (is.finite(grid_lo) && is.finite(grid_hi) && grid_lo <= grid_hi) {
    grid <- seq(grid_lo, grid_hi, length.out = 401)
    pred <- predict(best$fit, newdata = data.frame(x = grid, h = pmax(0, grid - best$cp)))
    idx <- which.min(pred)
    min_t <- grid[[idx]]
    min_y <- pred[[idx]]
  } else {
    min_t <- NA_real_
    min_y <- NA_real_
  }

  tibble(
    n_obs = n,
    n_unique_time = length(unique(round(d$x, 8))),
    change_point_weeks = best$cp,
    nearest_observed_to_change_weeks = nearest_obs,
    change_point_in_delivery_gap = in_gap,
    rss_linear = rss_linear,
    rss_piecewise = best$rss,
    delta_aic_vs_linear = aic_linear - best$aic,
    slope_before_change = slope_before,
    slope_after_change = slope_after,
    model_min_time_within_close_window = min_t,
    model_min_sws_within_close_window = min_y
  )
}

theme_ptb <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", color = ILLINI_BLUE),
      plot.subtitle = element_text(color = "grey30")
    )
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
data_dir <- file.path(project_root, "data")
out_dir <- Sys.getenv("PTB_SWS_OUT_DIR",
                      unset = file.path(project_root, "results", "SWS_delivery_changepoint_bland_altman"))
table_dir <- file.path(out_dir, "tables")
plot_dir <- file.path(out_dir, "plots")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

detected_k <- available_k_values(data_dir)
if (length(detected_k) == 0) stop("No V1_V7_longformat_k*.csv files found in ", data_dir)
k_values <- parse_k_values(detected_k)
k_values <- intersect(k_values, detected_k)
if (length(k_values) == 0) stop("Requested k values do not match available data files.")
reference_k <- as.integer(Sys.getenv("PTB_REFERENCE_K", as.character(max(k_values))))
if (!reference_k %in% k_values) {
  stop("PTB_REFERENCE_K=", reference_k, " is not in the selected k values.")
}
reference_k_value <- reference_k

message("Project root: ", project_root)
message("Data dir: ", data_dir)
message("Output dir: ", out_dir)
message("k values: ", paste(k_values, collapse = ", "), "; reference k=", reference_k)

dat <- bind_rows(lapply(k_values, read_long_k, data_dir = data_dir)) %>%
  mutate(
    Outcome_label = factor(Outcome_label, levels = c("FTB", "sPTB")),
    period = factor(as.character(period), levels = c("Prepartum", "Delivery", "Postpartum"))
  ) %>%
  arrange(k, Participant_ID, Visit_ID)

dat_all_groups <- bind_rows(
  dat %>% mutate(Outcome_group = as.character(Outcome_label)),
  dat %>% mutate(Outcome_group = "All")
) %>%
  mutate(Outcome_group = factor(Outcome_group, levels = c("All", "FTB", "sPTB")))

timing_by_visit <- dat %>%
  filter(k == min(k_values)) %>%
  group_by(Visit_ID, period) %>%
  summarise(
    n_obs = n(),
    n_subjects = n_distinct(Participant_ID),
    time_to_delivery_min_weeks = min(time_to_delivery_weeks),
    time_to_delivery_max_weeks = max(time_to_delivery_weeks),
    time_to_delivery_range_weeks = fmt_range(time_to_delivery_weeks),
    sws_min = min(SWS),
    sws_median = median(SWS),
    sws_max = max(SWS),
    .groups = "drop"
  )
write_csv(timing_by_visit, file.path(table_dir, "timing_by_visit_k_min.csv"))

delivery_gap_by_k <- dat %>%
  group_by(k) %>%
  summarise(
    n_obs = n(),
    closest_pre_delivery_weeks = max_or_na(time_to_delivery_weeks[time_to_delivery_weeks < 0]),
    closest_post_delivery_weeks = min_or_na(time_to_delivery_weeks[time_to_delivery_weeks > 0]),
    minimum_abs_time_to_delivery_weeks = min(abs(time_to_delivery_weeks)),
    has_delivery_measurement = any(abs(time_to_delivery_weeks) < 1e-8),
    .groups = "drop"
  )
write_csv(delivery_gap_by_k, file.path(table_dir, "delivery_gap_by_k.csv"))

close_dat <- dat_all_groups %>%
  filter(abs(time_to_delivery_weeks) <= CLOSE_WINDOW_WEEKS)

near_delivery_summary <- close_dat %>%
  group_by(k, Outcome_group) %>%
  summarise(
    n_obs = n(),
    n_subjects = n_distinct(Participant_ID),
    time_to_delivery_min_weeks = min(time_to_delivery_weeks),
    time_to_delivery_max_weeks = max(time_to_delivery_weeks),
    closest_pre_delivery_weeks = max_or_na(time_to_delivery_weeks[time_to_delivery_weeks < 0]),
    closest_post_delivery_weeks = min_or_na(time_to_delivery_weeks[time_to_delivery_weeks > 0]),
    sws_min = min(SWS),
    sws_q10 = q_or_na(SWS, 0.10),
    sws_q25 = q_or_na(SWS, 0.25),
    sws_median = median(SWS),
    sws_mean = mean(SWS),
    sws_q75 = q_or_na(SWS, 0.75),
    sws_q90 = q_or_na(SWS, 0.90),
    sws_max = max(SWS),
    .groups = "drop"
  )

observed_minimum_rows <- close_dat %>%
  group_by(k, Outcome_group) %>%
  filter(SWS == min(SWS, na.rm = TRUE)) %>%
  summarise(
    observed_min_sws = first(SWS),
    observed_min_time_min_weeks = min(time_to_delivery_weeks),
    observed_min_time_max_weeks = max(time_to_delivery_weeks),
    observed_min_visit_ids = paste(sort(unique(Visit_ID)), collapse = ";"),
    n_min_rows = n(),
    example_participant_ids = paste(head(sort(unique(Participant_ID)), 10), collapse = ";"),
    .groups = "drop"
  )

close_binned <- close_dat %>%
  add_time_bin(TIME_BIN_WIDTH_WEEKS) %>%
  group_by(k, Outcome_group, period, time_bin_left, time_bin_right, time_bin_mid, time_bin) %>%
  summarise(
    n_obs = n(),
    n_subjects = n_distinct(Participant_ID),
    sws_min = min(SWS),
    sws_q25 = q_or_na(SWS, 0.25),
    sws_median = median(SWS),
    sws_mean = mean(SWS),
    sws_q75 = q_or_na(SWS, 0.75),
    sws_max = max(SWS),
    .groups = "drop"
  )
write_csv(close_binned, file.path(table_dir, "near_delivery_sws_binned_summary.csv"))

minimum_median_bin <- close_binned %>%
  group_by(k, Outcome_group) %>%
  slice_min(order_by = sws_median, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(
    minimum_median_bin_left_weeks = time_bin_left,
    minimum_median_bin_right_weeks = time_bin_right,
    minimum_median_bin_mid_weeks = time_bin_mid,
    minimum_median_bin = time_bin,
    minimum_bin_n_obs = n_obs,
    minimum_bin_n_subjects = n_subjects,
    minimum_bin_sws_min = sws_min,
    minimum_bin_sws_q25 = sws_q25,
    minimum_bin_sws_median = sws_median,
    minimum_bin_sws_mean = sws_mean,
    minimum_bin_sws_q75 = sws_q75,
    minimum_bin_sws_max = sws_max
  )

near_delivery_minimum_range <- near_delivery_summary %>%
  left_join(observed_minimum_rows, by = c("k", "Outcome_group")) %>%
  left_join(minimum_median_bin, by = c("k", "Outcome_group"))
write_csv(near_delivery_minimum_range,
          file.path(table_dir, "near_delivery_minimum_sws_range_by_k.csv"))

change_input <- dat_all_groups %>%
  filter(abs(time_to_delivery_weeks) <= CHANGE_WINDOW_WEEKS)

change_point_results <- change_input %>%
  group_by(k, Outcome_group) %>%
  group_modify(~ fit_piecewise_change(
    .x,
    close_window = CLOSE_WINDOW_WEEKS,
    min_n = CHANGE_MIN_N,
    step = CHANGE_GRID_STEP_WEEKS
  )) %>%
  ungroup()
write_csv(change_point_results, file.path(table_dir, "sws_change_point_by_k.csv"))

trajectory_binned <- dat %>%
  add_time_bin(TIME_BIN_WIDTH_WEEKS) %>%
  group_by(k, Outcome_label, period, time_bin_left, time_bin_right, time_bin_mid, time_bin) %>%
  summarise(
    n_obs = n(),
    n_subjects = n_distinct(Participant_ID),
    sws_q25 = q_or_na(SWS, 0.25),
    sws_median = median(SWS),
    sws_mean = mean(SWS),
    sws_q75 = q_or_na(SWS, 0.75),
    .groups = "drop"
  )
write_csv(trajectory_binned, file.path(table_dir, "sws_time_to_delivery_binned_summary_by_k.csv"))

sptb_visit_summary <- dat %>%
  filter(Outcome_label == "sPTB") %>%
  group_by(k, Visit_ID, period) %>%
  summarise(
    n_obs = n(),
    n_subjects = n_distinct(Participant_ID),
    time_to_delivery_min_weeks = min(time_to_delivery_weeks),
    time_to_delivery_max_weeks = max(time_to_delivery_weeks),
    sws_q25 = q_or_na(SWS, 0.25),
    sws_median = median(SWS),
    sws_mean = mean(SWS),
    sws_q75 = q_or_na(SWS, 0.75),
    .groups = "drop"
  )
write_csv(sptb_visit_summary, file.path(table_dir, "sptb_visit_summary_by_k.csv"))

sptb_binned <- trajectory_binned %>%
  filter(Outcome_label == "sPTB")
write_csv(sptb_binned, file.path(table_dir, "sptb_time_to_delivery_binned_summary_by_k.csv"))

ref_dat <- dat %>%
  filter(k == reference_k) %>%
  select(
    Participant_ID,
    Visit_ID,
    Outcome_label,
    period,
    time_to_delivery_weeks,
    SWS_reference = SWS
  )

ba_reference <- bind_rows(lapply(setdiff(k_values, reference_k), function(k_comp) {
  comp <- dat %>%
    filter(k == k_comp) %>%
    select(Participant_ID, Visit_ID, SWS_comparison = SWS)

  ref_dat %>%
    inner_join(comp, by = c("Participant_ID", "Visit_ID")) %>%
    mutate(
      comparison_k = k_comp,
      reference_k = reference_k,
      comparison = paste0("k", comparison_k, " - k", reference_k),
      mean_sws = 0.5 * (SWS_comparison + SWS_reference),
      diff_sws = SWS_comparison - SWS_reference,
      abs_diff_sws = abs(diff_sws)
    )
}))

if (nrow(ba_reference) > 0) {
  write_csv(ba_reference, file.path(table_dir, "bland_altman_sws_k_vs_reference_pairs.csv"))

  ba_summary_period <- ba_reference %>%
    group_by(comparison_k, reference_k, comparison, period) %>%
    summarise(
      n_pairs = n(),
      bias = mean(diff_sws),
      sd_diff = sd(diff_sws),
      loa_lower = bias - 1.96 * sd_diff,
      loa_upper = bias + 1.96 * sd_diff,
      median_abs_diff = median(abs_diff_sws),
      q95_abs_diff = q_or_na(abs_diff_sws, 0.95),
      .groups = "drop"
    )
  write_csv(ba_summary_period,
            file.path(table_dir, "bland_altman_sws_k_vs_reference_by_period.csv"))

  ba_summary_outcome <- ba_reference %>%
    group_by(comparison_k, reference_k, comparison, period, Outcome_label) %>%
    summarise(
      n_pairs = n(),
      bias = mean(diff_sws),
      sd_diff = sd(diff_sws),
      loa_lower = bias - 1.96 * sd_diff,
      loa_upper = bias + 1.96 * sd_diff,
      median_abs_diff = median(abs_diff_sws),
      q95_abs_diff = q_or_na(abs_diff_sws, 0.95),
      .groups = "drop"
    )
  write_csv(ba_summary_outcome,
            file.path(table_dir, "bland_altman_sws_k_vs_reference_by_period_outcome.csv"))
}

pairwise_summary <- bind_rows(lapply(seq_along(k_values), function(i) {
  bind_rows(lapply(seq_along(k_values), function(j) {
    if (j <= i) return(NULL)
    k_a <- k_values[[i]]
    k_b <- k_values[[j]]
    a <- dat %>%
      filter(k == k_a) %>%
      select(Participant_ID, Visit_ID, Outcome_label, period, SWS_a = SWS)
    b <- dat %>%
      filter(k == k_b) %>%
      select(Participant_ID, Visit_ID, SWS_b = SWS)

    a %>%
      inner_join(b, by = c("Participant_ID", "Visit_ID")) %>%
      mutate(
        k_a = k_a,
        k_b = k_b,
        diff_sws = SWS_a - SWS_b,
        abs_diff_sws = abs(diff_sws)
      ) %>%
      group_by(k_a, k_b, period) %>%
      summarise(
        n_pairs = n(),
        bias = mean(diff_sws),
        sd_diff = sd(diff_sws),
        loa_lower = bias - 1.96 * sd_diff,
        loa_upper = bias + 1.96 * sd_diff,
        median_abs_diff = median(abs_diff_sws),
        q95_abs_diff = q_or_na(abs_diff_sws, 0.95),
        .groups = "drop"
      )
  }))
}))
write_csv(pairwise_summary, file.path(table_dir, "bland_altman_sws_all_pairwise_k_summary.csv"))

if (nrow(minimum_median_bin) > 0) {
  p_min <- minimum_median_bin %>%
    filter(Outcome_group %in% c("All", "FTB", "sPTB")) %>%
    ggplot(aes(x = k, y = minimum_bin_sws_median, color = Outcome_group, group = Outcome_group)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.2) +
    geom_errorbar(aes(ymin = minimum_bin_sws_q25, ymax = minimum_bin_sws_q75),
                  width = 0.14, linewidth = 0.5) +
    scale_color_manual(values = c(All = "grey35", FTB = ILLINI_BLUE, sPTB = ILLINI_ORANGE)) +
    scale_x_continuous(breaks = k_values) +
    labs(
      title = "Near-delivery SWS minimum range by k",
      subtitle = paste0("Minimum median ", TIME_BIN_WIDTH_WEEKS,
                        "-week bin within +/-", CLOSE_WINDOW_WEEKS,
                        " weeks of delivery"),
      x = "k first acquisitions averaged for V1/V2",
      y = "SWS in minimum median bin",
      color = "Group"
    ) +
    theme_ptb()
  ggsave(file.path(plot_dir, "near_delivery_minimum_sws_range_by_k.png"),
         p_min, width = 9, height = 5.5, dpi = 300)
}

if (nrow(change_point_results) > 0) {
  p_cp <- change_point_results %>%
    filter(Outcome_group %in% c("All", "FTB", "sPTB"), is.finite(change_point_weeks)) %>%
    ggplot(aes(x = k, y = change_point_weeks, color = Outcome_group, group = Outcome_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_line(linewidth = 0.8) +
    geom_point(aes(shape = change_point_in_delivery_gap), size = 2.4) +
    scale_color_manual(values = c(All = "grey35", FTB = ILLINI_BLUE, sPTB = ILLINI_ORANGE)) +
    scale_x_continuous(breaks = k_values) +
    labs(
      title = "Descriptive SWS change point by k",
      subtitle = paste0("Piecewise-linear fit within +/-", CHANGE_WINDOW_WEEKS,
                        " weeks of delivery; dashed line marks delivery"),
      x = "k first acquisitions averaged for V1/V2",
      y = "Estimated change point, weeks from delivery",
      color = "Group",
      shape = "In delivery gap"
    ) +
    theme_ptb()
  ggsave(file.path(plot_dir, "sws_change_point_by_k.png"),
         p_cp, width = 9, height = 5.5, dpi = 300)
}

p_traj <- trajectory_binned %>%
  filter(n_obs >= 2) %>%
  ggplot(aes(x = time_bin_mid, y = sws_median, color = Outcome_label, group = Outcome_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ k, ncol = 5) +
  scale_color_manual(values = c(FTB = ILLINI_BLUE, sPTB = ILLINI_ORANGE), drop = FALSE) +
  labs(
    title = "SWS trajectories on the time-to-delivery clock",
    subtitle = paste0("Binned medians; bin width = ", TIME_BIN_WIDTH_WEEKS,
                      " weeks; delivery is shown as x=0 but is not observed"),
    x = "Weeks from delivery",
    y = "Median SWS",
    color = "Outcome"
  ) +
  theme_ptb(base_size = 10)
ggsave(file.path(plot_dir, "sws_time_to_delivery_trajectory_by_k.png"),
       p_traj, width = 13, height = 7.5, dpi = 300)

if (nrow(sptb_binned) > 0) {
  k_palette <- grDevices::colorRampPalette(c(ILLINI_BLUE, ILLINI_ORANGE))(length(k_values))
  names(k_palette) <- as.character(k_values)
  p_sptb <- sptb_binned %>%
    filter(n_obs >= 2) %>%
    ggplot(aes(x = time_bin_mid, y = sws_median, color = factor(k), group = k)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.4) +
    scale_color_manual(values = k_palette) +
    labs(
      title = "sPTB SWS trajectories after time-to-delivery normalization",
      subtitle = paste0("Binned medians across k values; bin width = ",
                        TIME_BIN_WIDTH_WEEKS, " weeks"),
      x = "Weeks from delivery",
      y = "Median SWS among sPTB participants",
      color = "k"
    ) +
    theme_ptb()
  ggsave(file.path(plot_dir, "sptb_time_to_delivery_trajectory_by_k.png"),
         p_sptb, width = 9.5, height = 5.8, dpi = 300)
}

if (exists("ba_reference") && nrow(ba_reference) > 0) {
  ba_plot_dat <- ba_reference %>%
    mutate(comparison = factor(comparison,
                               levels = paste0("k", setdiff(k_values, reference_k_value),
                                               " - k", reference_k_value)))

  p_ba <- ggplot(ba_plot_dat, aes(x = mean_sws, y = diff_sws, color = Outcome_label)) +
    geom_hline(
      data = ba_summary_period,
      aes(yintercept = bias),
      inherit.aes = FALSE,
      color = ILLINI_BLUE,
      linewidth = 0.6
    ) +
    geom_hline(
      data = ba_summary_period,
      aes(yintercept = loa_lower),
      inherit.aes = FALSE,
      color = ILLINI_ORANGE,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    geom_hline(
      data = ba_summary_period,
      aes(yintercept = loa_upper),
      inherit.aes = FALSE,
      color = ILLINI_ORANGE,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    geom_point(alpha = 0.28, size = 0.75) +
    facet_grid(period ~ comparison, scales = "free_x") +
    scale_color_manual(values = c(FTB = ILLINI_BLUE, sPTB = ILLINI_ORANGE), drop = FALSE) +
    labs(
      title = paste0("Bland-Altman plots for SWS: k vs k", reference_k_value),
      subtitle = "Blue line is mean bias; orange dashed lines are +/-1.96 SD limits of agreement",
      x = "Mean SWS",
      y = paste0("SWS(k) - SWS(k", reference_k_value, ")"),
      color = "Outcome"
    ) +
    theme_ptb(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(plot_dir, "bland_altman_sws_k_vs_reference.png"),
         p_ba, width = 18, height = 7.5, dpi = 300)
}

if (nrow(pairwise_summary) > 0) {
  p_heat <- pairwise_summary %>%
    filter(period %in% c("Prepartum", "Postpartum")) %>%
    mutate(
      k_a = factor(k_a, levels = k_values),
      k_b = factor(k_b, levels = k_values)
    ) %>%
    ggplot(aes(x = k_a, y = k_b, fill = bias)) +
    geom_tile(color = "white", linewidth = 0.3) +
    facet_wrap(~ period) +
    scale_fill_gradient2(low = ILLINI_BLUE, mid = "white", high = ILLINI_ORANGE,
                         midpoint = 0) +
    labs(
      title = "Pairwise SWS bias across k values",
      subtitle = "Bias is SWS(k on x-axis) minus SWS(k on y-axis)",
      x = "k_a",
      y = "k_b",
      fill = "Bias"
    ) +
    theme_ptb()
  ggsave(file.path(plot_dir, "pairwise_k_sws_bias_heatmap.png"),
         p_heat, width = 8.5, height = 4.8, dpi = 300)
}

notes <- c(
  "SWS delivery-normalized analyses",
  "",
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Project root: ", project_root),
  paste0("Input files: ", paste(sprintf("V1_V7_longformat_k%d.csv", k_values), collapse = ", ")),
  paste0("Output directory: ", out_dir),
  "",
  "Clock definition:",
  "  time_to_delivery_weeks = GA_weeks - GA_delivery_weeks",
  "  negative values are prepartum; positive values are postpartum.",
  "  A delivery scan is not imputed. The observed delivery gap is reported in tables/delivery_gap_by_k.csv.",
  "",
  "Near-delivery minimum:",
  paste0("  close window = +/-", CLOSE_WINDOW_WEEKS, " weeks from delivery."),
  paste0("  binned minimum uses ", TIME_BIN_WIDTH_WEEKS, "-week bins and selects the bin with lowest median SWS."),
  "",
  "Change point:",
  paste0("  descriptive piecewise-linear change points are fit within +/-",
         CHANGE_WINDOW_WEEKS, " weeks from delivery."),
  "  change_point_in_delivery_gap marks estimates falling between the closest pre- and postpartum observations.",
  "",
  "Bland-Altman:",
  paste0("  reference k = ", reference_k, "."),
  "  The current preprocessing changes V1/V2 with k. Postpartum V3-V7 values are averaged over all acquisitions,",
  "  so postpartum cross-k Bland-Altman differences should be zero unless preprocessing is changed.",
  "",
  "Plots are PNG only; no PDF plot output is generated."
)
writeLines(notes, con = file.path(out_dir, "analysis_notes.txt"))

manifest <- tibble(
  artifact = c(
    "analysis_notes.txt",
    file.path("tables", "timing_by_visit_k_min.csv"),
    file.path("tables", "delivery_gap_by_k.csv"),
    file.path("tables", "near_delivery_minimum_sws_range_by_k.csv"),
    file.path("tables", "near_delivery_sws_binned_summary.csv"),
    file.path("tables", "sws_change_point_by_k.csv"),
    file.path("tables", "sws_time_to_delivery_binned_summary_by_k.csv"),
    file.path("tables", "sptb_visit_summary_by_k.csv"),
    file.path("tables", "sptb_time_to_delivery_binned_summary_by_k.csv"),
    file.path("tables", "bland_altman_sws_k_vs_reference_pairs.csv"),
    file.path("tables", "bland_altman_sws_k_vs_reference_by_period.csv"),
    file.path("tables", "bland_altman_sws_k_vs_reference_by_period_outcome.csv"),
    file.path("tables", "bland_altman_sws_all_pairwise_k_summary.csv"),
    file.path("plots", "near_delivery_minimum_sws_range_by_k.png"),
    file.path("plots", "sws_change_point_by_k.png"),
    file.path("plots", "sws_time_to_delivery_trajectory_by_k.png"),
    file.path("plots", "sptb_time_to_delivery_trajectory_by_k.png"),
    file.path("plots", "bland_altman_sws_k_vs_reference.png"),
    file.path("plots", "pairwise_k_sws_bias_heatmap.png")
  )
) %>%
  mutate(
    path = file.path(out_dir, artifact),
    exists = file.exists(path)
  )
write_csv(manifest, file.path(out_dir, "manifest.csv"))

message("Done. Wrote outputs to: ", out_dir)
