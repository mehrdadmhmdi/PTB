#!/usr/bin/env Rscript

# FTB-only semiparametric phase-variance estimator for all biomarkers and taus.
# This is a new analysis script.  It does not estimate individual b_i values.
# It estimates sigma_b^2 using within-subject quantile-hit covariances.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(quantreg)
  library(splines)
  library(tidyr)
  library(knitr)
  library(parallel)
})

DATA_PATH <- Sys.getenv("PTB_DATA_PATH", unset = "data/V1_V7_longformat_k1.csv")
OUT_DIR <- Sys.getenv("PTB_SEMIPARAM_OUT_DIR", unset = "Test_results_semiparametric_FTB")
B_BOOT <- as.integer(Sys.getenv("B_BOOT", unset = "200"))
if (is.na(B_BOOT) || B_BOOT < 0) B_BOOT <- 200
NCORES <- as.integer(Sys.getenv("NCORES", unset = "1"))
if (is.na(NCORES) || NCORES < 1) NCORES <- 1

VARS <- c("SWS", "Midband", "Intercept", "AC", "Slope")
TAUS <- c(0.20, 0.50, 0.80)
METHODS <- c("quadratic", "spline")
DF_SPLINE <- 3
H_SPARSE <- 0.05
SEED0 <- 2026

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

prep_ftb <- function(path) {
  raw <- read_csv(path, show_col_types = FALSE)
  ga_raw <- suppressWarnings(as.numeric(raw$GA))
  ga_mult <- if (max(ga_raw, na.rm = TRUE) <= 2) 39 else 1

  raw %>%
    mutate(
      GA = as.numeric(GA) * ga_mult,
      Outcome = as.character(Outcome),
      Participant_ID = as.character(Participant_ID)
    ) %>%
    filter(Outcome == "0") %>%
    arrange(Participant_ID, GA)
}

fit_quad <- function(dat, yvar, tau) {
  rq(as.formula(paste0(yvar, " ~ GA + I(GA^2)")), tau = tau, data = dat)
}

pred_quad <- function(fit, t) {
  b <- coef(fit)
  as.numeric(b[1] + b[2] * t + b[3] * t^2)
}

deriv_quad <- function(fit, t) {
  b <- coef(fit)
  as.numeric(b[2] + 2 * b[3] * t)
}

fit_spline <- function(dat, yvar, tau, df_spline = 3) {
  f <- as.formula(paste0(yvar, " ~ ns(GA, df = ", df_spline, ")"))
  rq(f, tau = tau, data = dat)
}

pred_rq <- function(fit, t) {
  as.numeric(predict(fit, newdata = data.frame(GA = t)))
}

deriv_spline <- function(fit, t, delta = 0.05) {
  (pred_rq(fit, t + delta) - pred_rq(fit, t - delta)) / (2 * delta)
}

pair_contrib <- function(d, lambda_tol = 1e-8) {
  if (nrow(d) < 2) return(data.frame(N = 0, D = 0, npair = 0))

  cmb <- combn(seq_len(nrow(d)), 2)
  ll <- d$lambda[cmb[1, ]] * d$lambda[cmb[2, ]]
  hh <- d$H[cmb[1, ]] * d$H[cmb[2, ]]
  keep <- is.finite(ll) & is.finite(hh) & abs(ll) > lambda_tol

  data.frame(
    N = sum(ll[keep] * hh[keep]),
    D = sum(ll[keep]^2),
    npair = sum(keep)
  )
}

estimate_sigma_b <- function(dat, yvar, tau = 0.5, h = 0.05,
                             method = c("quadratic", "spline"),
                             df_spline = 3, min_sparsity = 1e-4) {
  method <- match.arg(method)

  dat <- dat %>%
    filter(
      !is.na(.data[[yvar]]),
      !is.na(GA),
      !is.na(Participant_ID)
    )

  if (nrow(dat) == 0) stop("No rows after filtering for ", yvar)

  if (method == "quadratic") {
    fit0 <- fit_quad(dat, yvar, tau)
    fitm <- fit_quad(dat, yvar, tau - h)
    fitp <- fit_quad(dat, yvar, tau + h)
    q0 <- pred_quad(fit0, dat$GA)
    qm <- pred_quad(fitm, dat$GA)
    qp <- pred_quad(fitp, dat$GA)
    dq <- deriv_quad(fit0, dat$GA)
  } else {
    fit0 <- fit_spline(dat, yvar, tau, df_spline)
    fitm <- fit_spline(dat, yvar, tau - h, df_spline)
    fitp <- fit_spline(dat, yvar, tau + h, df_spline)
    q0 <- pred_rq(fit0, dat$GA)
    qm <- pred_rq(fitm, dat$GA)
    qp <- pred_rq(fitp, dat$GA)
    dq <- deriv_spline(fit0, dat$GA)
  }

  sparsity <- pmax((qp - qm) / (2 * h), min_sparsity)
  ghat <- 1 / sparsity

  work <- dat %>%
    mutate(
      H = as.numeric(.data[[yvar]] <= q0) - tau,
      lambda = ghat * dq
    ) %>%
    group_by(Participant_ID) %>%
    group_modify(~ pair_contrib(.x)) %>%
    ungroup()

  sig2 <- sum(work$N) / sum(work$D)

  list(
    sigma_b2 = sig2,
    sigma_b2_positive = max(0, sig2),
    n_pairs = sum(work$npair),
    n = nrow(dat),
    n_id = length(unique(dat$Participant_ID))
  )
}

bootstrap_sigma_b <- function(dat, yvar, tau, method, B = 200,
                              h = 0.05, df_spline = 3, seed = 1) {
  ids <- unique(dat$Participant_ID)
  n_ids <- length(ids)
  rows_by_id <- split(seq_len(nrow(dat)), dat$Participant_ID)

  set.seed(seed)
  boot <- replicate(B, {
    sid <- sample(ids, size = n_ids, replace = TRUE)
    db <- bind_rows(lapply(seq_along(sid), function(m) {
      d <- dat[rows_by_id[[sid[m]]], , drop = FALSE]
      d$Participant_ID <- paste0(sid[m], "__boot__", m)
      d
    }))

    out <- try(
      estimate_sigma_b(
        db,
        yvar = yvar,
        tau = tau,
        h = h,
        method = method,
        df_spline = df_spline
      )$sigma_b2,
      silent = TRUE
    )

    if (inherits(out, "try-error")) NA_real_ else out
  })

  boot[is.finite(boot)]
}

one_run <- function(dat, yvar, tau, method, B) {
  point <- estimate_sigma_b(
    dat,
    yvar = yvar,
    tau = tau,
    h = H_SPARSE,
    method = method,
    df_spline = DF_SPLINE
  )

  if (B > 0) {
    boot <- bootstrap_sigma_b(
      dat,
      yvar = yvar,
      tau = tau,
      method = method,
      B = B,
      h = H_SPARSE,
      df_spline = DF_SPLINE,
      seed = SEED0 + as.integer(100 * tau) + match(yvar, VARS) * 10 + match(method, METHODS)
    )
    ci <- as.numeric(quantile(boot, c(0.025, 0.975), na.rm = TRUE))
    se <- sd(boot)
    p_norm <- if (is.finite(se) && se > 0) 2 * pnorm(-abs(point$sigma_b2 / se)) else NA_real_
    p_upper_boot <- mean(boot <= 0)
    b_ok <- length(boot)
  } else {
    ci <- c(NA_real_, NA_real_)
    se <- NA_real_
    p_norm <- NA_real_
    p_upper_boot <- NA_real_
    b_ok <- 0
  }

  data.frame(
    biomarker = yvar,
    tau = tau,
    method = method,
    sigma_b2 = point$sigma_b2,
    sigma_b2_positive = point$sigma_b2_positive,
    boot_se = se,
    ci_low = ci[1],
    ci_high = ci[2],
    p_two_sided_normal = p_norm,
    p_upper_boot = p_upper_boot,
    B_ok = b_ok,
    n_pairs = point$n_pairs,
    n = point$n,
    n_id = point$n_id,
    stringsAsFactors = FALSE
  )
}

dat_ftb <- prep_ftb(DATA_PATH)
message("FTB rows: ", nrow(dat_ftb))
message("FTB subjects: ", length(unique(dat_ftb$Participant_ID)))
message("Bootstrap replicates per row: ", B_BOOT)
message("Parallel workers across biomarker-tau-method rows: ", NCORES)

grid <- expand.grid(
  biomarker = VARS,
  tau = TAUS,
  method = METHODS,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

run_grid_row <- function(k) {
  message(
    "Running ", k, "/", nrow(grid),
    ": ", grid$biomarker[k],
    " tau=", grid$tau[k],
    " method=", grid$method[k]
  )
  one_run(dat_ftb, grid$biomarker[k], grid$tau[k], grid$method[k], B_BOOT)
}

if (NCORES > 1) {
  cl <- makeCluster(NCORES)
  on.exit(stopCluster(cl), add = TRUE)

  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(readr)
      library(dplyr)
      library(quantreg)
      library(splines)
      library(tidyr)
    })
    NULL
  })

  clusterExport(
    cl,
    varlist = c(
      "grid", "dat_ftb", "B_BOOT", "H_SPARSE", "DF_SPLINE", "SEED0",
      "VARS", "TAUS", "METHODS",
      "fit_quad", "pred_quad", "deriv_quad",
      "fit_spline", "pred_rq", "deriv_spline",
      "pair_contrib", "estimate_sigma_b", "bootstrap_sigma_b",
      "one_run", "run_grid_row"
    ),
    envir = environment()
  )

  res <- bind_rows(parLapply(cl, seq_len(nrow(grid)), run_grid_row))
} else {
  res <- bind_rows(lapply(seq_len(nrow(grid)), run_grid_row))
}

out_csv <- file.path(OUT_DIR, "ftb_semiparametric_phase_variance_all.csv")
write_csv(res, out_csv)

fmt_tex_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(
    is.na(x),
    "--",
    ifelse(abs(x) < 0.0005, "0.000", sprintf("%.3f", x))
  )
}

res_tex <- res %>%
  mutate(
    sigma_b2 = fmt_tex_num(sigma_b2),
    sigma_b2_positive = fmt_tex_num(sigma_b2_positive),
    ci = ifelse(
      is.na(ci_low),
      "--",
      paste0("(", fmt_tex_num(ci_low), ", ", fmt_tex_num(ci_high), ")")
    ),
    p_two_sided_normal = fmt_tex_num(p_two_sided_normal),
    p_upper_boot = fmt_tex_num(p_upper_boot)
  ) %>%
  select(biomarker, tau, method, sigma_b2, sigma_b2_positive, ci, p_two_sided_normal, p_upper_boot, B_ok)

out_tex <- file.path(OUT_DIR, "ftb_semiparametric_phase_variance_all.tex")
sink(out_tex)
cat("\\begin{table*}[!htbp]\n")
cat("\\centering\n")
cat("\\caption{FTB-only semiparametric phase-variance estimates for all biomarkers. All rows use 939 observations from 454 FTB subjects. Each row uses subject-level bootstrap resamples.}\n")
cat("\\label{tab:ftb_semiparametric_phase_variance_summary}\n")
cat("\\scriptsize\n")
print(
  knitr::kable(
    res_tex,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    col.names = c("Biomarker", "$\\tau$", "Method", "$\\widehat\\sigma_b^2$", "$\\widehat\\sigma_{b,+}^2$", "95\\% CI", "Normal $p$", "$p_{\\le 0}^{\\mathrm{boot}}$", "$B_{ok}$")
  )
)
cat("\\end{table*}\n")
sink()

message("Wrote ", out_csv)
message("Wrote ", out_tex)
