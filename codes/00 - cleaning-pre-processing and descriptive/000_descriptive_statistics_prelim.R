rm(list = ls())
start_time <- Sys.time()

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(tidyr)
})

options(digits = 3)

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    if (!is.null(active_path) && nzchar(active_path)) {
      return(dirname(normalizePath(active_path, mustWork = FALSE)))
    }
  }

  getwd()
}

mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
sd_na <- function(x) if (sum(!is.na(x)) < 2) NA_real_ else sd(x, na.rm = TRUE)
median_na <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
q_na <- function(x, p) if (all(is.na(x))) NA_real_ else unname(quantile(x, p, na.rm = TRUE))
min_na <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
max_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

first_existing <- function(data, vars) {
  intersect(vars, names(data))
}

is_yes <- function(x) {
  value <- tolower(trimws(as.character(x)))
  value %in% c("1", "yes", "y", "true")
}

is_no <- function(x) {
  value <- tolower(trimws(as.character(x)))
  value %in% c("0", "no", "n", "false")
}

scaled_to_weeks <- function(x) {
  if (all(is.na(x))) {
    return(x)
  }

  if (max(x, na.rm = TRUE) <= 5) {
    return(x * 39)
  }

  x
}

summarise_continuous <- function(data, vars, group_vars = character()) {
  vars <- first_existing(data, vars)
  if (length(vars) == 0) {
    return(tibble())
  }

  data %>%
    pivot_longer(all_of(vars), names_to = "variable", values_to = "value") %>%
    group_by(across(all_of(group_vars)), variable) %>%
    summarise(
      N = sum(!is.na(value)),
      Missing = sum(is.na(value)),
      Mean = mean_na(value),
      SD = sd_na(value),
      Median = median_na(value),
      Q1 = q_na(value, 0.25),
      Q3 = q_na(value, 0.75),
      Min = min_na(value),
      Max = max_na(value),
      .groups = "drop"
    )
}

summarise_categorical <- function(data, vars, group_vars = character()) {
  vars <- first_existing(data, vars)
  if (length(vars) == 0) {
    return(tibble())
  }

  data %>%
    pivot_longer(all_of(vars), names_to = "variable", values_to = "level",
                 values_transform = list(level = as.character)) %>%
    mutate(level = if_else(is.na(level) | level == "", "Missing", level)) %>%
    count(across(all_of(group_vars)), variable, level, name = "N") %>%
    group_by(across(all_of(group_vars)), variable) %>%
    mutate(
      Denominator = sum(N),
      Percent = 100 * N / Denominator
    ) %>%
    ungroup()
}

format_iqr <- function(median, q1, q3) {
  ifelse(
    is.na(median),
    NA_character_,
    sprintf("%.1f [%.1f, %.1f]", median, q1, q3)
  )
}

theme_ptb <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
}

script_dir <- get_script_dir()
data_path <- Sys.getenv(
  "PTB_DATA_PATH",
  unset = file.path(script_dir, "..", "..", "data", "V1_V7_longformat_k1.csv")
)
data_path <- normalizePath(data_path, mustWork = TRUE)

output_dir <- Sys.getenv(
  "PTB_DESC_OUTPUT_DIR",
  unset = file.path(script_dir, "descriptive_statistics")
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
plot_dir <- file.path(output_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(plot, filename, width = 8, height = 5) {
  path <- file.path(plot_dir, filename)
  ggsave(path, plot, width = width, height = height, dpi = 300, bg = "white")
  path
}

raw <- read_csv(data_path, show_col_types = FALSE)

df <- raw %>%
  mutate(
    GA_weeks = scaled_to_weeks(GA),
    GA_delivery_weeks = scaled_to_weeks(GA_delivery),
    is_mptb = is_yes(mPTB),
    ptb_status = case_when(
      is_yes(Outcome) ~ "Preterm",
      is_no(Outcome) ~ "Term",
      !is.na(GA_delivery_weeks) & GA_delivery_weeks < 37 ~ "Preterm",
      !is.na(GA_delivery_weeks) & GA_delivery_weeks >= 37 ~ "Term",
      TRUE ~ NA_character_
    ),
    ptb_status = factor(ptb_status, levels = c("Term", "Preterm")),
    Visit_ID = as.integer(Visit_ID)
  ) %>%
  filter(!is_mptb)

participant <- df %>%
  arrange(Participant_ID, Visit_ID) %>%
  group_by(Participant_ID) %>%
  summarise(
    across(
      any_of(c(
        "Age", "race", "Ethnicity", "pregnancy_count",
        "last_pregnancy_outcome", "FTB_count", "prior_PTB",
        "prior_PTB_count", "prior_abortion", "abortion_count",
        "smoke_cigarettes", "drink_alcohol", "drugs_used",
        "GA_delivery_weeks", "Outcome", "mPTB", "baby_sex",
        "blood_pressure", "diabetes", "liver_disease", "lung_disease",
        "high_cholesterol", "anemia", "ptb_status"
      )),
      first
    ),
    n_visits = n_distinct(Visit_ID),
    first_visit_week = min_na(GA_weeks),
    last_visit_week = max_na(GA_weeks),
    .groups = "drop"
  )

visit_counts <- df %>%
  count(Visit_ID, ptb_status, name = "observations") %>%
  group_by(Visit_ID) %>%
  mutate(percent_within_visit = 100 * observations / sum(observations)) %>%
  ungroup()

visits_per_participant <- participant %>%
  count(n_visits, name = "participants") %>%
  mutate(percent = 100 * participants / sum(participants))

overview <- tibble(
  Metric = c(
    "Input file",
    "Raw observations",
    "Analyzed observations after excluding medically indicated PTB",
    "Raw participants",
    "Analyzed participants",
    "Term participants",
    "Spontaneous preterm participants",
    "Medically indicated PTB participants excluded",
    "Median visits per analyzed participant",
    "Maximum visits per analyzed participant"
  ),
  Value = c(
    data_path,
    as.character(nrow(raw)),
    as.character(nrow(df)),
    as.character(n_distinct(raw$Participant_ID)),
    as.character(n_distinct(df$Participant_ID)),
    as.character(sum(participant$ptb_status == "Term", na.rm = TRUE)),
    as.character(sum(participant$ptb_status == "Preterm", na.rm = TRUE)),
    as.character(n_distinct(raw$Participant_ID[is_yes(raw$mPTB)])),
    as.character(median(participant$n_visits, na.rm = TRUE)),
    as.character(max(participant$n_visits, na.rm = TRUE))
  )
)

participant_continuous_vars <- c(
  "Age", "pregnancy_count", "FTB_count", "prior_PTB_count",
  "abortion_count", "GA_delivery_weeks", "n_visits",
  "first_visit_week", "last_visit_week"
)

participant_categorical_vars <- c(
  "race", "Ethnicity", "last_pregnancy_outcome", "prior_PTB",
  "prior_abortion", "smoke_cigarettes", "drink_alcohol", "drugs_used",
  "baby_sex", "blood_pressure", "diabetes", "liver_disease",
  "lung_disease", "high_cholesterol", "anemia"
)

visit_continuous_vars <- c(
  "GA_weeks", "CL", "AC", "Slope", "Intercept", "Midband",
  "Kappa", "Mu", "SWS", "ESDG", "EACG", "ESDS", "EACS"
)

participant_continuous <- summarise_continuous(
  participant,
  participant_continuous_vars,
  group_vars = "ptb_status"
)

participant_continuous_overall <- summarise_continuous(
  participant,
  participant_continuous_vars
) %>%
  mutate(ptb_status = "Overall", .before = 1)

participant_categorical <- summarise_categorical(
  participant,
  participant_categorical_vars,
  group_vars = "ptb_status"
)

participant_categorical_overall <- summarise_categorical(
  participant,
  participant_categorical_vars
) %>%
  mutate(ptb_status = "Overall", .before = 1)

biomarker_by_visit <- summarise_continuous(
  df,
  visit_continuous_vars,
  group_vars = c("Visit_ID", "ptb_status")
)

biomarker_by_visit_overall <- summarise_continuous(
  df,
  visit_continuous_vars,
  group_vars = "Visit_ID"
) %>%
  mutate(ptb_status = "Overall", .after = "Visit_ID")

missingness <- df %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "missing_n") %>%
  mutate(
    total_n = nrow(df),
    missing_percent = 100 * missing_n / total_n
  ) %>%
  arrange(desc(missing_percent), variable)

ga_visit_summary <- biomarker_by_visit_overall %>%
  filter(variable == "GA_weeks") %>%
  transmute(
    Visit_ID,
    N,
    Median_GA_weeks = Median,
    GA_IQR_weeks = format_iqr(Median, Q1, Q3),
    Min_GA_weeks = Min,
    Max_GA_weeks = Max
  )

plot_files <- character()
add_plot <- function(plot, filename, width = 8, height = 5) {
  plot_files <<- c(plot_files, save_plot(plot, filename, width, height))
}

ILLINI_ORANGE <- "#FF5F05"
ILLINI_BLUE <- "#13294B"
outcome_colors <- c("Term" = ILLINI_BLUE, "Preterm" = ILLINI_ORANGE)

participant_plot <- participant %>%
  filter(!is.na(ptb_status))

visit_plot <- df %>%
  filter(!is.na(ptb_status)) %>%
  mutate(Visit = factor(paste0("V", Visit_ID),
                        levels = paste0("V", sort(unique(Visit_ID)))))

if (nrow(participant_plot) > 0) {
  p_visits_per_participant <- ggplot(
    participant_plot,
    aes(x = factor(n_visits), fill = ptb_status)
  ) +
    geom_bar(position = "stack", color = "white", linewidth = 0.2) +
    scale_fill_manual(values = outcome_colors, drop = FALSE) +
    labs(
      title = "Number of Observed Visits Per Participant",
      x = "Observed visits",
      y = "Participants",
      fill = "Outcome"
    ) +
    theme_ptb()

  add_plot(p_visits_per_participant, "plot_01_visits_per_participant.png")
}

if (nrow(visit_plot) > 0) {
  p_visit_counts <- ggplot(
    visit_plot,
    aes(x = Visit, fill = ptb_status)
  ) +
    geom_bar(color = "white", linewidth = 0.2) +
    scale_fill_manual(values = outcome_colors, drop = FALSE) +
    labs(
      title = "Observed Records By Visit",
      x = "Visit",
      y = "Observations",
      fill = "Outcome"
    ) +
    theme_ptb()

  add_plot(p_visit_counts, "plot_02_observations_by_visit.png")

  p_ga_visit <- ggplot(
    visit_plot,
    aes(x = Visit, y = GA_weeks, fill = ptb_status)
  ) +
    geom_boxplot(outlier.alpha = 0.25, width = 0.65, position = position_dodge(width = 0.75)) +
    scale_fill_manual(values = outcome_colors, drop = FALSE) +
    labs(
      title = "Gestational Age At Each Visit",
      x = "Visit",
      y = "Gestational age at scan (weeks)",
      fill = "Outcome"
    ) +
    theme_ptb()

  add_plot(p_ga_visit, "plot_03_ga_by_visit.png", width = 8, height = 5)
}

if (nrow(participant_plot) > 0) {
  p_ga_delivery <- ggplot(
    participant_plot,
    aes(x = GA_delivery_weeks, fill = ptb_status)
  ) +
    geom_histogram(alpha = 0.75, bins = 24, position = "identity",
                   color = "white") +
    geom_vline(xintercept = 37, color = ILLINI_BLUE,
               linetype = "dashed", linewidth = 0.5) +
    facet_wrap(~ ptb_status, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = outcome_colors, drop = FALSE) +
    labs(
      title = "Gestational Age At Delivery",
      subtitle = "Dashed line marks 37 weeks",
      x = "Gestational age at delivery (weeks)",
      y = "Participants",
      fill = "Outcome"
    ) +
    theme_ptb()

  add_plot(p_ga_delivery, "plot_04_ga_delivery_by_outcome.png", width = 8, height = 6)
}

if (nrow(visit_plot) > 0) {
  p_ga_continuous <- ggplot(
    visit_plot %>% filter(!is.na(GA_weeks)),
    aes(x = GA_weeks, fill = ptb_status, color = ptb_status)
  ) +
    geom_density(alpha = 0.22, linewidth = 0.7, na.rm = TRUE) +
    scale_fill_manual(values = outcome_colors, drop = FALSE) +
    scale_color_manual(values = outcome_colors, drop = FALSE) +
    labs(
      title = "Continuous Gestational Age Distribution At Scans",
      x = "Gestational age at scan (weeks)",
      y = "Density",
      fill = "Outcome",
      color = "Outcome"
    ) +
    theme_ptb()

  add_plot(p_ga_continuous, "plot_08_scan_ga_density_continuous.png",
           width = 8, height = 5)
}

biomarker_plot_vars <- first_existing(
  df,
  c("CL", "AC", "SWS", "Kappa", "Mu", "Midband", "Slope", "Intercept")
)

if (length(biomarker_plot_vars) > 0 && nrow(visit_plot) > 0) {
  biomarker_long <- visit_plot %>%
    select(Visit_ID, Visit, GA_weeks, ptb_status, all_of(biomarker_plot_vars)) %>%
    pivot_longer(all_of(biomarker_plot_vars),
                 names_to = "variable", values_to = "value") %>%
    filter(!is.na(value), !is.na(ptb_status))

  if (nrow(biomarker_long) > 0) {
    p_biomarker_box <- ggplot(
      biomarker_long,
      aes(x = Visit, y = value, fill = ptb_status)
    ) +
      geom_boxplot(outlier.alpha = 0.2, width = 0.65, position = position_dodge(width = 0.75)) +
      facet_wrap(~ variable, scales = "free_y", ncol = 2) +
      scale_fill_manual(values = outcome_colors, drop = FALSE) +
      labs(
        title = "QUS And Cervical-Length Distributions By Visit",
        x = "Visit",
        y = "Observed value",
        fill = "Outcome"
      ) +
      theme_ptb()

    add_plot(p_biomarker_box, "plot_05_biomarker_distributions_by_visit.png",
             width = 10, height = 9)

    biomarker_trend <- biomarker_long %>%
      group_by(Visit_ID, Visit, ptb_status, variable) %>%
      summarise(
        Median = median_na(value),
        Q1 = q_na(value, 0.25),
        Q3 = q_na(value, 0.75),
        .groups = "drop"
      )

    p_biomarker_trend <- ggplot(
      biomarker_trend,
      aes(x = Visit_ID, y = Median, color = ptb_status, group = ptb_status)
    ) +
      geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.08, alpha = 0.65) +
      geom_line(linewidth = 0.55) +
      geom_point(size = 1.8) +
      facet_wrap(~ variable, scales = "free_y", ncol = 2) +
      scale_x_continuous(breaks = sort(unique(biomarker_trend$Visit_ID)),
                         labels = paste0("V", sort(unique(biomarker_trend$Visit_ID)))) +
      scale_color_manual(values = outcome_colors, drop = FALSE) +
      labs(
        title = "Median Biomarker Trajectories",
        subtitle = "Points show medians; bars show interquartile ranges",
        x = "Visit",
        y = "Median value",
        color = "Outcome"
      ) +
      theme_ptb()

    add_plot(p_biomarker_trend, "plot_06_biomarker_median_iqr_trends.png",
             width = 10, height = 9)

    continuous_ga_vars <- first_existing(
      df,
      c("CL", "AC", "SWS", "Kappa", "Mu", "Midband")
    )

    biomarker_ga <- biomarker_long %>%
      filter(variable %in% continuous_ga_vars, !is.na(GA_weeks))

    if (nrow(biomarker_ga) > 0) {
      p_biomarker_ga <- ggplot(
        biomarker_ga,
        aes(x = GA_weeks, y = value, color = ptb_status)
      ) +
        geom_point(alpha = 0.28, size = 0.8) +
        geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                    linewidth = 0.75, na.rm = TRUE) +
        facet_wrap(~ variable, scales = "free_y", ncol = 2) +
        scale_color_manual(values = outcome_colors, drop = FALSE) +
        labs(
          title = "Biomarkers Across Continuous Gestational Age",
          subtitle = "Curves are loess smooths over scan gestational age",
          x = "Gestational age at scan (weeks)",
          y = "Observed value",
          color = "Outcome"
        ) +
        theme_ptb()

      add_plot(p_biomarker_ga, "plot_09_biomarker_vs_ga_continuous.png",
               width = 10, height = 8)
    }
  }
}

top_missing <- missingness %>%
  filter(missing_n > 0) %>%
  slice_max(order_by = missing_percent, n = 15, with_ties = FALSE) %>%
  mutate(variable = reorder(variable, missing_percent))

if (nrow(top_missing) > 0) {
  p_missing <- ggplot(
    top_missing,
    aes(x = variable, y = missing_percent)
  ) +
    geom_col(fill = ILLINI_BLUE, width = 0.7) +
    coord_flip() +
    scale_y_continuous(labels = label_percent(scale = 1)) +
    labs(
      title = "Variables With Highest Missingness",
      x = NULL,
      y = "Missing observations"
    ) +
    theme_ptb() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

  add_plot(p_missing, "plot_07_missingness.png", width = 8, height = 5.5)
}

write_csv(overview, file.path(output_dir, "00_cohort_overview.csv"))
write_csv(visits_per_participant, file.path(output_dir, "01_visits_per_participant.csv"))
write_csv(visit_counts, file.path(output_dir, "02_visit_counts_by_outcome.csv"))
write_csv(
  bind_rows(participant_continuous_overall, participant_continuous),
  file.path(output_dir, "03_participant_continuous_by_outcome.csv")
)
write_csv(
  bind_rows(participant_categorical_overall, participant_categorical),
  file.path(output_dir, "04_participant_categorical_by_outcome.csv")
)
write_csv(
  bind_rows(biomarker_by_visit_overall, biomarker_by_visit),
  file.path(output_dir, "05_visit_biomarker_summary_by_outcome.csv")
)
write_csv(ga_visit_summary, file.path(output_dir, "06_ga_by_visit.csv"))
write_csv(missingness, file.path(output_dir, "07_missingness.csv"))

summary_file <- file.path(output_dir, "descriptive_summary.txt")
sink(summary_file)
cat("PTB V1-V7 Descriptive Summary\n")
cat("=============================\n\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input:", data_path, "\n")
cat("Output directory:", output_dir, "\n\n")

cat("Cohort Overview\n")
cat("---------------\n")
print(overview, n = nrow(overview))

cat("\nVisits Per Participant\n")
cat("----------------------\n")
print(visits_per_participant, n = nrow(visits_per_participant))

cat("\nVisit Counts By Outcome\n")
cat("-----------------------\n")
print(visit_counts, n = nrow(visit_counts))

cat("\nGestational Age By Visit\n")
cat("------------------------\n")
print(ga_visit_summary, n = nrow(ga_visit_summary))

cat("\nParticipant-Level Continuous Variables By Outcome\n")
cat("-------------------------------------------------\n")
print(participant_continuous, n = nrow(participant_continuous))

cat("\nVisit-Level Biomarkers By Visit And Outcome\n")
cat("-------------------------------------------\n")
print(biomarker_by_visit, n = nrow(biomarker_by_visit))

cat("\nHighest Missingness\n")
cat("-------------------\n")
print(head(missingness, 25), n = 25)

cat("\nGenerated Plots\n")
cat("---------------\n")
if (length(plot_files) == 0) {
  cat("No plots were generated.\n")
} else {
  cat(paste0("  - ", basename(plot_files)), sep = "\n")
  cat("\n")
}
sink()

elapsed_time <- round(difftime(Sys.time(), start_time, units = "secs"), 2)
cat("Descriptive statistics complete in", elapsed_time, "seconds.\n")
cat("Input:", data_path, "\n")
cat("Analyzed observations:", nrow(df), "\n")
cat("Analyzed participants:", n_distinct(df$Participant_ID), "\n")
cat("Output directory:", output_dir, "\n")
cat("Main report:", summary_file, "\n")
cat("Plot directory:", plot_dir, "\n")
