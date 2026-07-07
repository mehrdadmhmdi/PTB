#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Quantile regression with natural splines in GA
# FTB-only analysis
#
# Model fitted for each biomarker v and tau:
#   Q_tau(Y_ij^(v) | GA_ij) = alpha_tau^(v) + B(GA_ij)^T theta_tau^(v)
#
# Test:
#   H0: theta_tau^(v) = 0
#   i.e. no GA effect in the FTB group; the quantile trajectory is constant over GA
#
# Outputs:
#   Test_results_GA_FTB/qr_models_FTB_baseOnly_GA.rds
#   Test_results_GA_FTB/qr_models_FTB_baseOnly_coef_long_GA.csv
#   Test_results_GA_FTB/qr_models_FTB_baseOnly_coef_wide_GA.csv
#   Test_results_GA_FTB/qr_FTB_constantOverTime_clusterBoot_GA.csv
#   Test_results_GA_FTB/qr_FTB_constantOverTime_clusterBoot_GA.tex
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

# Keep only FTB observations
df$Outcome <- as.factor(df$Outcome)

if ("0" %in% levels(df$Outcome)) {
  df <- df %>% filter(Outcome == "0")
} else {
  stop("FTB was not found among Outcome levels.")
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

dir.create("Test_results_GA_FTB", showWarnings = FALSE, recursive = TRUE)

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

coef_to_df <- function(mod, var, tau) {
  cf <- coef(mod)
  data.frame(
    var      = var,
    tau      = tau,
    term     = names(cf),
    estimate = as.numeric(cf),
    n_rows   = get_nobs_rq(mod),
    obj      = qr_obj(mod, tau),
    stringsAsFactors = FALSE
  )
}

make_formula <- function(v, df_spline) {
  as.formula(sprintf("%s ~ ns(GA, df=%d)", v, df_spline))
}

get_ns_main_terms <- function(coef_names, df_spline) {
  pat_ns_main <- sprintf("^ns\\(GA\\s*,\\s*df\\s*=\\s*%d\\)", df_spline)
  coef_names[grepl(pat_ns_main, coef_names)]
}

safe_inverse <- function(M) {
  tryCatch(solve(M), error = function(e) qr.solve(M))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3) Fit and store FTB-only base models
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MODEL_LIST <- list()
COEF_LIST  <- list()

for (v in vars2) {
  for (t in taus) {
    f_base <- make_formula(v, df_spline)
    m_base <- rq(f_base, tau = t, data = df)
    
    key <- sprintf("%s__tau_%s", v, format(t, scientific = FALSE, trim = TRUE))
    MODEL_LIST[[key]] <- m_base
    
    COEF_LIST[[length(COEF_LIST) + 1]] <- coef_to_df(m_base, v, t)
  }
}

coef_long <- bind_rows(COEF_LIST) %>%
  arrange(var, tau, term)

write.csv(
  coef_long,
  "Test_results_GA_FTB/qr_models_FTB_baseOnly_coef_long_GA.csv",
  row.names = FALSE
)

coef_wide <- coef_long %>%
  select(var, tau, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  arrange(var, tau)

write.csv(
  coef_wide,
  "Test_results_GA_FTB/qr_models_FTB_baseOnly_coef_wide_GA.csv",
  row.names = FALSE
)

saveRDS(
  MODEL_LIST,
  "Test_results_GA_FTB/qr_models_FTB_baseOnly_GA.rds"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4) Subject-level cluster bootstrap Wald test for:
#    H0: theta = 0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
idx_by_id <- split(seq_len(nrow(df)), df[[ID_COL]])
uniq_ids  <- names(idx_by_id)
n_ids     <- length(uniq_ids)

ftb_constant_over_time_test_cluster_boot <- function(v, tau, B = 500, seed = 1) {
  f_base <- make_formula(v, df_spline)
  
  fit0 <- rq(f_base, tau = tau, data = df)
  cn   <- names(coef(fit0))
  
  test_names <- get_ns_main_terms(cn, df_spline)
  
  if (length(test_names) == 0) {
    return(data.frame(
      var       = v,
      tau       = tau,
      df_test   = 0,
      W_obs     = NA_real_,
      p_boot    = NA_real_,
      p_chisq   = NA_real_,
      n_boot_ok = 0,
      note      = "No spline terms found.",
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
    
    fb <- try(rq(f_base, tau = tau, data = dat_b), silent = TRUE)
    
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

out_ftb_constant_time <- bind_rows(
  lapply(vars2, function(v) {
    bind_rows(
      lapply(taus, function(t) {
        ftb_constant_over_time_test_cluster_boot(
          v    = v,
          tau  = t,
          B    = B_boot,
          seed = SEED0 + 4000 + round(1000 * t) + match(v, vars2)
        )
      })
    )
  })
) %>%
  relocate(var, tau)

write.csv(
  out_ftb_constant_time,
  "Test_results_GA_FTB/qr_FTB_constantOverTime_clusterBoot_GA.csv",
  row.names = FALSE
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5) LaTeX table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tab_ftb_constant_time <- out_ftb_constant_time %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  arrange(var, tau)

latex_ftb_constant_time <- kable(
  tab_ftb_constant_time,
  format   = "latex",
  booktabs = TRUE,
  caption  = paste0(
    "FTB-only cluster-bootstrap Wald test for constant over time. ",
    "$H_0$: all GA spline coefficients are zero, so the biomarker quantile ",
    "trajectory is constant over gestational age among FTB pregnancies."
  ),
  label    = "tab:qr_ftb_constant_over_time_clusterboot",
  escape   = FALSE
)

writeLines(
  latex_ftb_constant_time,
  "Test_results_GA_FTB/qr_FTB_constantOverTime_clusterBoot_GA.tex"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6) Optional: print LaTeX table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n\n===== LaTeX table: FTB-only constant over time =====\n")
cat(latex_ftb_constant_time)
cat("\n")
