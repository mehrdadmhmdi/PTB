#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Quantile regression natural-spline models (base + full) for multiple biomarkers and taus
# + Subject-level (cluster) bootstrap Wald test for GLOBAL time-dependence:
#
#   H0^{Global time}: (theta_{tau1},...,theta_{tauM}, gamma_{tau1},...,gamma_{tauM}) = 0
#   i.e., no spline/time pattern in either group (FTB or sPTB); allows a constant level shift delta_tau.
#
# Outputs:
#   Test_results/qr_models_all_coef_long.csv
#   Test_results/qr_models_all_coef_wide.csv
#   Test_results/qr_models_all.rds
#   Test_results/qr_globalTime_thetaPlusGamma_clusterBoot.csv
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())
library(readr)
library(dplyr)
library(tidyr)
library(quantreg)
library(splines)
library(knitr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data loading + preprocessing (as in your script)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- read_csv("../../data/V1_V7_longformat_k1.csv")

df <- data
df <- df %>% mutate(GA = GA * 39, GA_delivery = GA_delivery * 39)

# Make Outcome factor with explicit reference = FTB (adjust levels if your labels differ)
df$Outcome <- as.factor(df$Outcome)
if (all(c("FTB","sPTB") %in% levels(df$Outcome))) {
  df$Outcome <- factor(df$Outcome, levels = c("FTB","sPTB"))
} else {
  # If your labels differ, keep as-is, but you should set reference explicitly.
  # Example: df$Outcome <- relevel(df$Outcome, ref = "FTB_label_here")
  message("Outcome levels are: ", paste(levels(df$Outcome), collapse = ", "),
          ". Consider setting FTB as reference with relevel()/factor(levels=...).")
}

df$GA_diff <- df$GA - df$GA_delivery

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Settings
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vars2     <- c("SWS", "Midband", "Intercept", "AC", "Slope")
taus      <- c(0.2, 0.5, 0.8)
df_spline <- 3
B_boot    <- 500
SEED0     <- 2026

# Subject ID
ID_COL       <- "Participant_ID"
df[[ID_COL]] <- as.factor(df[[ID_COL]])

# output dir
dir.create("Test_results", showWarnings = FALSE, recursive = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Helper: objective value (check loss) for rq
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qr_obj <- function(mod, tau){
  u <- resid(mod)
  sum(u * (tau - (u < 0)))
}

# robust nobs for rq objects
get_nobs_rq <- function(mod){
  if (!is.null(mod$n))  return(as.integer(mod$n))
  if (!is.null(mod$X))  return(nrow(mod$X))
  if (!is.null(mod$x))  return(nrow(mod$x))
  r <- tryCatch(resid(mod), error = function(e) NULL)
  if (!is.null(r))      return(length(r))
  NA_integer_
}

coef_to_df <- function(mod, var, tau, model_name){
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Fit + store base/full models and export all coefficients
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MODEL_LIST <- list()
COEF_LIST  <- list()

for (v in vars2) {
  for (t in taus) {
    f_base <- as.formula(sprintf("%s ~ ns(GA_diff, df=%d)", v, df_spline))
    f_full <- as.formula(sprintf("%s ~ Outcome * ns(GA_diff, df=%d)", v, df_spline))
    
    m0 <- rq(f_base, tau = t, data = df)
    m1 <- rq(f_full, tau = t, data = df)
    
    key <- sprintf("%s__tau_%s", v, format(t, scientific = FALSE, trim = TRUE))
    MODEL_LIST[[key]] <- list(base = m0, full = m1)
    
    COEF_LIST[[length(COEF_LIST) + 1]] <- coef_to_df(m0, v, t, "base")
    COEF_LIST[[length(COEF_LIST) + 1]] <- coef_to_df(m1, v, t, "full")
  }
}

coef_long <- bind_rows(COEF_LIST) %>%
  arrange(var, tau, model, term)

write.csv(coef_long, "Test_results/qr_models_all_coef_long.csv", row.names = FALSE)
saveRDS(MODEL_LIST, "Test_results/qr_models_all.rds")

coef_wide <- coef_long %>%
  select(var, tau, model, term, estimate) %>%
  pivot_wider(names_from = model, values_from = estimate) %>%
  mutate(delta_full_minus_base = full - base) %>%
  arrange(var, tau, term)

write.csv(coef_wide, "Test_results/qr_models_all_coef_wide.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2) GLOBAL time-dependence cluster-bootstrap Wald test:
#    H0: theta = 0 AND gamma = 0  (allow delta != 0)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
idx_by_id <- split(seq_len(nrow(df)), df[[ID_COL]])
uniq_ids  <- names(idx_by_id)
n_ids     <- length(uniq_ids)

global_time_test_cluster_boot <- function(v, tau, B = 500, seed = 1) {
  f_full <- as.formula(sprintf("%s ~ Outcome * ns(GA_diff, df=%d)", v, df_spline))
  
  fit0 <- rq(f_full, tau = tau, data = df)
  cn   <- names(coef(fit0))
  
  # theta block: spline main effects (ns terms)
  # gamma block: Outcome:spline interactions
  pat_ns_main <- sprintf("^ns\\(GA_diff\\s*,\\s*df\\s*=\\s*%d\\)", df_spline)
  ns_main <- cn[grepl(pat_ns_main, cn)]
  ns_int  <- cn[grepl("^Outcome", cn) & grepl(":ns\\(GA_diff", cn)]
  
  test_names <- c(ns_main, ns_int)
  if (length(test_names) == 0) {
    return(data.frame(
      var=v, tau=tau, df_test=0, W_obs=NA_real_, p_boot=NA_real_, p_chisq=NA_real_,
      n_boot_ok=0, note="No spline terms found (check df_spline / formula / naming)."
    ))
  }
  
  beta_hat <- coef(fit0)[test_names]
  boot_mat <- matrix(NA_real_, nrow = B, ncol = length(test_names),
                     dimnames = list(NULL, test_names))
  
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
  
  boot_ctr <- sweep(boot_mat, 2, beta_hat, "-")
  V <- cov(boot_ctr, use = "pairwise.complete.obs") + diag(1e-10, length(beta_hat))
  
  W_obs  <- as.numeric(t(beta_hat) %*% solve(V, beta_hat))
  W_star <- apply(boot_ctr, 1, function(bv) {
    if (anyNA(bv)) return(NA_real_)
    as.numeric(t(bv) %*% solve(V, bv))
  })
  
  data.frame(
    var=v, tau=tau,
    df_test = length(beta_hat),
    W_obs = W_obs,
    p_boot = mean(W_star >= W_obs, na.rm = TRUE),
    p_chisq = pchisq(W_obs, df = length(beta_hat), lower.tail = FALSE),
    n_boot_ok = sum(is.finite(W_star)),
    note = ""
  )
}

out_global_time <- bind_rows(
  lapply(vars2, function(v) {
    bind_rows(lapply(taus, function(t) {
      global_time_test_cluster_boot(
        v, t, B = B_boot,
        seed = SEED0 + 4000 + round(1000*t) + match(v, vars2)
      )
    }))
  })
) %>% relocate(var, tau)

write.csv(out_global_time,
          "Test_results/qr_globalTime_thetaPlusGamma_clusterBoot.csv",
          row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (Optional) LaTeX table for the GLOBAL time-dependence test
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tab_time <- out_global_time %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  arrange(var, tau)

kable(tab_time, format = "latex", booktabs = TRUE,
      caption = "Cluster-bootstrap Wald test for global time-dependence: $H_0: \\theta=0$ and $\\gamma=0$ (natural spline df=3).",
      label = "tab:qr_global_time_clusterboot")
