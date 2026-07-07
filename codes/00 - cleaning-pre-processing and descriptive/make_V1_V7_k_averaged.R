## Build V1-V7 analysis files for multiple k values.
##
## Inputs, relative to git/data by default:
##   - V1_V7_wideformat_k1.csv: clinical/output scaffold used downstream
##   - raw/UIC_QUS_8.8.22_C425updated.xlsx: raw V1/V2 QUS acquisitions
##   - raw/postpartum_data.csv: raw V3-V7 postpartum QUS acquisitions
##
## Outputs:
##   - V1_V7_wideformat_k{k}.csv
##   - V1_V7_longformat_k{k}.csv
##
## Optional output:
##   - QUS_means_k.csv, if PTB_WRITE_QUS_MEANS=true

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(tidyr)
})

K_VALUES <- 1:10
QUS_SHEET <- "per-image (SLD)"
POSTPARTUM_AVERAGE_FIRST_K <- FALSE
WRITE_QUS_MEANS <- tolower(Sys.getenv("PTB_WRITE_QUS_MEANS", unset = "false")) %in%
  c("1", "true", "yes", "y")

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

parse_k_values <- function(default) {
  value <- Sys.getenv("PTB_K_VALUES", unset = "")
  if (!nzchar(value)) {
    return(default)
  }

  if (grepl("^\\s*\\d+\\s*:\\s*\\d+\\s*$", value)) {
    bounds <- as.integer(strsplit(gsub("\\s+", "", value), ":", fixed = TRUE)[[1]])
    return(seq(bounds[[1]], bounds[[2]]))
  }

  values <- as.integer(strsplit(gsub("\\s+", "", value), ",", fixed = TRUE)[[1]])
  values <- values[!is.na(values)]
  if (length(values) == 0) {
    stop("PTB_K_VALUES did not contain any valid integers.")
  }

  sort(unique(values))
}

mean_na <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

find_input <- function(filename, directories) {
  candidates <- file.path(directories, filename)
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop("Missing input file: ", filename, "\nChecked:\n  ",
         paste(candidates, collapse = "\n  "))
  }

  normalizePath(existing[[1]], mustWork = TRUE)
}

average_first_k <- function(data, k, group_cols, value_cols, keep_cols = character()) {
  data %>%
    group_by(across(all_of(group_cols))) %>%
    arrange(.acquisition_order, .row_order, .by_group = TRUE) %>%
    slice_head(n = k) %>%
    summarise(
      across(all_of(keep_cols), ~ first(.x)),
      across(all_of(value_cols), mean_na),
      .groups = "drop"
    )
}

script_dir <- get_script_dir()
data_dir <- normalizePath(file.path(script_dir, "..", "..", "data"), mustWork = TRUE)
raw_dir <- Sys.getenv("PTB_RAW_DIR", unset = file.path(data_dir, "raw"))
raw_dir <- normalizePath(raw_dir, mustWork = TRUE)
output_dir <- Sys.getenv("PTB_OUTPUT_DIR", unset = data_dir)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
k_values <- parse_k_values(K_VALUES)

template_path <- file.path(data_dir, "V1_V7_wideformat_k1.csv")
qus_path <- find_input("UIC_QUS_8.8.22_C425updated.xlsx", c(raw_dir, data_dir))
postpartum_path <- find_input("postpartum_data.csv", c(raw_dir, data_dir))

wide_template <- read_csv(template_path, show_col_types = FALSE)
qus_raw <- read_excel(qus_path, sheet = QUS_SHEET)
postpartum_raw <- read_csv(postpartum_path, show_col_types = FALSE, na = c("", "NULL"))

clinical_cols <- names(wide_template)[!grepl("_V\\d+$", names(wide_template))]
pre_visit_cols <- intersect(c("Participant_ID", "GA_V1", "GA_V2", "CL_V1", "CL_V2"),
                            names(wide_template))

clinical_scaffold <- wide_template %>%
  select(all_of(clinical_cols)) %>%
  left_join(wide_template %>% select(all_of(pre_visit_cols)), by = "Participant_ID")

pre_qus <- qus_raw %>%
  mutate(
    .row_order = row_number(),
    .acquisition_order = suppressWarnings(as.integer(Acquisition_index))
  ) %>%
  filter(Visit_ID %in% c("V1", "V2")) %>%
  rename(
    AC = AC_at_center_freq_in_dB_per_cm_MHz,
    Slope = BSC_LF_slope_in_dB_per_MHz,
    Intercept = BSC_LF_intercept_in_dB,
    Midband = BSC_LF_midband_fit_in_dB,
    ESDG = Gaussian_ESD_mean_in_um,
    EACG = Gaussian_EAC_mean_in_dB_per_cubic_cm,
    ESDS = Fluid_sphere_ESD_mean_in_um,
    EACS = Fluid_sphere_EAC_mean_in_dB_per_cubic_cm,
    Kappa = Envelope_stat_k_mean,
    Mu = Envelope_stat_mu_mean,
    SWS = SWS_mean
  )

pre_value_cols <- c("AC", "Slope", "Intercept", "Midband",
                    "Kappa", "Mu", "SWS", "ESDG", "EACG", "ESDS", "EACS")

postpartum <- postpartum_raw %>%
  mutate(
    .row_order = row_number(),
    .acquisition_order = suppressWarnings(as.integer(Acquisition_index)),
    GA_f_Del = Days_f_Del / (7 * 39)
  ) %>%
  rename(
    AC = AC_at_center_freq_in_dB_per_cm_MHz,
    Slope = BSC_LF_slope_in_dB_per_MHz,
    Intercept = BSC_LF_intercept_in_dB,
    Midband = `10Log10(BSC)`,
    SWS = SWS_mean
  )

post_value_cols <- c("GA_f_Del", "AC", "Slope", "Intercept", "Midband", "SWS")

for (k in k_values) {
  pre_avg <- average_first_k(
    data = pre_qus,
    k = k,
    group_cols = c("Participant_ID", "Visit_ID"),
    value_cols = pre_value_cols,
    keep_cols = c("Sonographer", "Repeat_ID")
  )

  qus_means <- pre_avg %>%
    transmute(
      Participant_ID,
      Visit_ID,
      Sonographer,
      Repeat_ID,
      AC,
      LF_Intercept = Intercept,
      LF_Midband = Midband,
      LF_Slope = Slope,
      Kappa,
      Mu,
      SWS,
      index = row_number()
    )

  if (WRITE_QUS_MEANS) {
    write_csv(qus_means, file.path(output_dir, sprintf("QUS_means_%d.csv", k)))
  }

  pre_wide <- pre_avg %>%
    select(Participant_ID, Visit_ID, all_of(pre_value_cols)) %>%
    pivot_wider(
      id_cols = Participant_ID,
      names_from = Visit_ID,
      values_from = all_of(pre_value_cols),
      names_glue = "{.value}_{Visit_ID}"
    )

  if (POSTPARTUM_AVERAGE_FIRST_K) {
    post_avg <- average_first_k(
      data = postpartum,
      k = k,
      group_cols = c("Participant_ID", "Visit_ID"),
      value_cols = post_value_cols
    )
  } else {
    post_avg <- postpartum %>%
      group_by(Participant_ID, Visit_ID) %>%
      summarise(across(all_of(post_value_cols), mean_na), .groups = "drop")
  }

  post_wide <- post_avg %>%
    pivot_wider(
      id_cols = Participant_ID,
      names_from = Visit_ID,
      values_from = all_of(post_value_cols),
      names_glue = "{.value}_{Visit_ID}"
    ) %>%
    rename_with(~ sub("^GA_f_Del_V", "GA_V", .x), starts_with("GA_f_Del_V"))

  df_wide <- clinical_scaffold %>%
    left_join(pre_wide, by = "Participant_ID") %>%
    left_join(post_wide, by = "Participant_ID") %>%
    mutate(across(starts_with("GA_V") & !matches("^GA_V[12]$"),
                  ~ .x + GA_delivery))

  wide_path <- file.path(output_dir, sprintf("V1_V7_wideformat_k%d.csv", k))
  write_csv(df_wide, wide_path)

  df_long <- df_wide %>%
    pivot_longer(
      cols = matches("_V\\d+$"),
      names_to = c(".value", "Visit_ID"),
      names_pattern = "^(.*)_V(\\d+)$"
    ) %>%
    mutate(Visit_ID = as.integer(Visit_ID)) %>%
    filter(!if_all(any_of(pre_value_cols), is.na))

  long_path <- file.path(output_dir, sprintf("V1_V7_longformat_k%d.csv", k))
  write_csv(df_long, long_path)

  cat(sprintf(
    "k=%d: wrote %s (%d x %d) and %s (%d x %d)\n",
    k,
    basename(wide_path), nrow(df_wide), ncol(df_wide),
    basename(long_path), nrow(df_long), ncol(df_long)
  ))
}
