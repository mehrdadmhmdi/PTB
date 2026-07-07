#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Quantile regression with natural splines in GA
# - Biomarkers: SWS, Midband, Intercept, AC, Slope
# - Quantiles: 0.2, 0.5, 0.8
#
# Models fitted for each biomarker and tau:
#   1) base     : y ~ ns(GA, df = df_spline)
#   2) additive : y ~ Outcome + ns(GA, df = df_spline)
#   3) full     : y ~ Outcome * ns(GA, df = df_spline)
#
# Wald tests with subject-level (cluster) bootstrap:
#   Test 1: H0 = trajectories are constant over time
#           theta = 0 and gamma = 0
#           i.e. no GA-spline effect in either group; a constant group shift is allowed
#
#   Test 2: H0 = no group effect in the GA trajectory
#           delta = 0 and gamma = 0
#           i.e. FTB and sPTB have the same entire quantile trajectory over GA
#
# Outputs:
#   Test_results_GA/qr_models_all_GA.rds
#   Test_results_GA/qr_models_all_coef_long_GA.csv
#   Test_results_GA/qr_models_all_coef_wide_GA.csv
#   Test_results_GA/qr_constantOverTime_clusterBoot_GA.csv
#   Test_results_GA/qr_noGroupEffect_clusterBoot_GA.csv
#   Test_results_GA/qr_constantOverTime_clusterBoot_GA.tex
#   Test_results_GA/qr_noGroupEffect_clusterBoot_GA.tex
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

library(readr)
library(dplyr)
library(tidyr)
library(quantreg)
library(splines)
library(knitr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 0) Data loading + preprocessing
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- read_csv("../../data/V1_V7_longformat_k1.csv", show_col_types = FALSE)

df <- data %>%
  mutate(
    GA = GA * 39,
    GA_delivery = GA_delivery * 39,
    GA_diff = GA - GA_delivery
  )

# Make Outcome factor with FTB as reference if available
df$Outcome <- as.factor(df$Outcome)
if (all(c("FTB", "sPTB") %in% levels(df$Outcome))) {
  df$Outcome <- factor(df$Outcome, levels = c("FTB", "sPTB"))
} else {
  message(
    "Outcome levels are: ", paste(levels(df$Outcome), collapse = ", "),
    ". If needed, manually set FTB as the reference."
  )
}

# Subject ID
ID_COL <- "Participant_ID"
df[[ID_COL]] <- as.factor(df[[ID_COL]])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Settings
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vars2     <- c("SWS", "Midband", "Intercept", "AC", "Slope")
taus      <- c(0.2, 0.5, 0.8)
df_spline <- 3
B_boot    <- 500
SEED0     <- 2026

dir.create("Test_results_GA", showWarnings = FALSE, recursive = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2) Helpers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qr_obj <- function(mod, tau) {
  u <- resid(mod)
  sum(u * (tau - (u < 0)))
}

get_nobs_rq <- function(mod) {
  if (!is.null(mod$n)) return(as.integer(mod$n))
  if (!is.null(mod$X)) return(nrow(mod$X))
  if (!is.null(mod$x)) return(nrow(mod$x))
  r <- tryCatch(resid(mod), error = function(e) NULL)
  if (!is.null(r)) return(length(r))
  NA_integer_
}

coef_to_df <- function(mod, var, tau, model_name) {
  cf <- coef(mod)
  data.frame(
    var      = var,
    tau      = tau,
    model    = model_name,
    term     = names(cf),
    estimate = as.numeric(cf),
    n_rows   = get_nobs_rq(mod),
    obj      = qr_obj(mod, tau),
    stringsAsFactors = FALSE
  )
}

make_formulas <- function(v, df_spline) {
  list(
    base     = as.formula(sprintf("%s ~ ns(GA, df=%d)", v, df_spline)),
    additive = as.formula(sprintf("%s ~ Outcome + ns(GA, df=%d)", v, df_spline)),
    full     = as.formula(sprintf("%s ~ Outcome * ns(GA, df=%d)", v, df_spline))
  )
}

get_ns_main_terms <- function(coef_names, df_spline) {
  pat_ns_main <- sprintf("^ns\\(GA\\s*,\\s*df\\s*=\\s*%d\\)", df_spline)
  coef_names[grepl(pat_ns_main, coef_names)]
}

get_ns_interaction_terms <- function(coef_names) {
  coef_names[grepl("^Outcome", coef_names) & grepl(":ns\\(GA", coef_names)]
}

get_outcome_main_terms <- function(coef_names) {
  coef_names[grepl("^Outcome", coef_names) & !grepl(":ns\\(GA", coef_names)]
}

safe_inverse <- function(M) {
  tryCatch(solve(M), error = function(e) qr.solve(M))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3) Fit and store models
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MODEL_LIST <- list()
COEF_LIST  <- list()

for (v in vars2) {
  for (t in taus) {
    forms <- make_formulas(v, df_spline)
    
    m_base <- rq(forms$base, tau = t, data = df)
    m_add  <- rq(forms$additive, tau = t, data = df)
    m_full <- rq(forms$full, tau = t, data = df)
    
    key <- sprintf("%s__tau_%s", v, format(t, scientific = FALSE, trim = TRUE))
    MODEL_LIST[[key]] <- list(
      base     = m_base,
      additive = m_add,
      full     = m_full
    )
    
    COEF_LIST[[length(COEF_LIST) + 1]] <- coef_to_df(m_base, v, t, "base")
    COEF_LIST[[length(COEF_LIST) + 1]] <- coef_to_df(m_add,  v, t, "additive")
    COEF_LIST[[length(COEF_LIST) + 1]] <- coef_to_df(m_full, v, t, "full")
  }
}

coef_long <- bind_rows(COEF_LIST) %>%
  arrange(var, tau, model, term)

write.csv(
  coef_long,
  "Test_results_GA/qr_models_all_coef_long_GA.csv",
  row.names = FALSE
)

saveRDS(
  MODEL_LIST,
  "Test_results_GA/qr_models_all_GA.rds"
)

coef_wide <- coef_long %>%
  select(var, tau, model, term, estimate) %>%
  pivot_wider(names_from = model, values_from = estimate) %>%
  mutate(
    delta_additive_minus_base = additive - base,
    delta_full_minus_additive = full - additive,
    delta_full_minus_base     = full - base
  ) %>%
  arrange(var, tau, term)

write.csv(
  coef_wide,
  "Test_results_GA/qr_models_all_coef_wide_GA.csv",
  row.names = FALSE
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4) Cluster bootstrap setup
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
idx_by_id <- split(seq_len(nrow(df)), df[[ID_COL]])
uniq_ids  <- names(idx_by_id)
n_ids     <- length(uniq_ids)

block_test_cluster_boot <- function(v, tau, term_selector, B = 500, seed = 1, note_label = "") {
  f_full <- as.formula(sprintf("%s ~ Outcome * ns(GA, df=%d)", v, df_spline))
  
  fit0 <- rq(f_full, tau = tau, data = df)
  cn   <- names(coef(fit0))
  
  test_names <- term_selector(cn)
  
  if (length(test_names) == 0) {
    return(data.frame(
      var       = v,
      tau       = tau,
      df_test   = 0,
      W_obs     = NA_real_,
      p_boot    = NA_real_,
      p_chisq   = NA_real_,
      n_boot_ok = 0,
      note      = paste("No terms found for", note_label),
      stringsAsFactors = FALSE
    ))
  }
  
  beta_hat <- coef(fit0)[test_names]
  
  boot_mat <- matrix(
    NA_real_,
    nrow = B,
    ncol = length(test_names),
    dimnames = list(NULL, test_names)
  )
  
  set.seed(seed)
  
  for (b in seq_len(B)) {
    samp_ids <- sample(uniq_ids, size = n_ids, replace = TRUE)
    rows_b   <- unlist(idx_by_id[samp_ids], use.names = FALSE)
    dat_b    <- df[rows_b, , drop = FALSE]
    
    fb <- try(rq(f_full, tau = tau, data = dat_b), silent = TRUE)
    
    if (!inherits(fb, "try-error")) {
      cb <- coef(fb)
      boot_mat[b, ] <- cb[test_names]
    }
  }
  
  ok <- complete.cases(boot_mat)
  n_boot_ok <- sum(ok)
  
  if (n_boot_ok < 10) {
    return(data.frame(
      var       = v,
      tau       = tau,
      df_test   = length(test_names),
      W_obs     = NA_real_,
      p_boot    = NA_real_,
      p_chisq   = NA_real_,
      n_boot_ok = n_boot_ok,
      note      = "Too few successful bootstrap replicates.",
      stringsAsFactors = FALSE
    ))
  }
  
  boot_ctr <- sweep(boot_mat[ok, , drop = FALSE], 2, beta_hat, "-")
  V <- cov(boot_ctr)
  
  if (is.null(dim(V))) {
    V <- matrix(V, nrow = 1, ncol = 1)
  }
  
  V <- V + diag(1e-10, ncol(V))
  Vinv <- safe_inverse(V)
  
  W_obs <- as.numeric(t(beta_hat) %*% Vinv %*% beta_hat)
  
  W_star <- apply(boot_ctr, 1, function(bv) {
    as.numeric(t(bv) %*% Vinv %*% bv)
  })
  
  data.frame(
    var       = v,
    tau       = tau,
    df_test   = length(beta_hat),
    W_obs     = W_obs,
    p_boot    = mean(W_star >= W_obs, na.rm = TRUE),
    p_chisq   = pchisq(W_obs, df = length(beta_hat), lower.tail = FALSE),
    n_boot_ok = n_boot_ok,
    note      = "",
    stringsAsFactors = FALSE
  )
}

# Test 1: H0: theta = 0 and gamma = 0
constant_term_selector <- function(cn) {
  c(
    get_ns_main_terms(cn, df_spline),
    get_ns_interaction_terms(cn)
  )
}

# Test 2: H0: delta = 0 and gamma = 0
group_effect_term_selector <- function(cn) {
  c(
    get_outcome_main_terms(cn),
    get_ns_interaction_terms(cn)
  )
}

run_test_grid <- function(term_selector, seed_offset, note_label) {
  bind_rows(
    lapply(vars2, function(v) {
      bind_rows(
        lapply(taus, function(t) {
          block_test_cluster_boot(
            v             = v,
            tau           = t,
            term_selector = term_selector,
            B             = B_boot,
            seed          = SEED0 + seed_offset + round(1000 * t) + match(v, vars2),
            note_label    = note_label
          )
        })
      )
    })
  ) %>%
    relocate(var, tau)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5) Test 1: H0 = constant over time
#    theta = 0 and gamma = 0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out_constant_time <- run_test_grid(
  term_selector = constant_term_selector,
  seed_offset   = 4000,
  note_label    = "constant-over-time test"
)

write.csv(
  out_constant_time,
  "Test_results_GA/qr_constantOverTime_clusterBoot_GA.csv",
  row.names = FALSE
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6) Test 2: H0 = no group effect in the GA trajectory
#    delta = 0 and gamma = 0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out_group_effect <- run_test_grid(
  term_selector = group_effect_term_selector,
  seed_offset   = 8000,
  note_label    = "no-group-effect test"
)

write.csv(
  out_group_effect,
  "Test_results_GA/qr_noGroupEffect_clusterBoot_GA.csv",
  row.names = FALSE
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7) LaTeX tables for both tests
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tab_constant_time <- out_constant_time %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  arrange(var, tau)

tab_group_effect <- out_group_effect %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  arrange(var, tau)

latex_constant_time <- kable(
  tab_constant_time,
  format   = "latex",
  booktabs = TRUE,
  caption  = paste0(
    "Cluster-bootstrap Wald test for constant over time. ",
    "$H_0$: all GA spline main effects and all Outcome$\\times$GA spline ",
    "interaction effects are zero; both groups are constant over GA, while ",
    "a constant group shift is allowed."
  ),
  label    = "tab:qr_constant_over_time_clusterboot",
  escape   = FALSE
)

latex_group_effect <- kable(
  tab_group_effect,
  format   = "latex",
  booktabs = TRUE,
  caption  = paste0(
    "Cluster-bootstrap Wald test for no group effect in the GA trajectory. ",
    "$H_0$: the Outcome main effect and all Outcome$\\times$GA spline ",
    "interaction effects are zero; FTB and sPTB have the same entire ",
    "quantile trajectory over GA."
  ),
  label    = "tab:qr_no_group_effect_clusterboot",
  escape   = FALSE
)

writeLines(
  latex_constant_time,
  "Test_results_GA/qr_constantOverTime_clusterBoot_GA.tex"
)

writeLines(
  latex_group_effect,
  "Test_results_GA/qr_noGroupEffect_clusterBoot_GA.tex"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 8) Optional: print tables to console
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n\n===== LaTeX table: constant over time =====\n")
cat(latex_constant_time)
cat("\n\n===== LaTeX table: no group effect =====\n")
cat(latex_group_effect)
cat("\n")