#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(quantreg)
  library(splines)
})

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
    n_id = length(unique(dat$Participant_ID)),
    fit = fit0
  )
}

bootstrap_sigma_b <- function(dat, yvar, tau = 0.5, h = 0.05,
                              method = c("quadratic", "spline"),
                              df_spline = 3, B = 200, seed = 2026) {
  method <- match.arg(method)
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

  boot <- boot[is.finite(boot)]
  point <- estimate_sigma_b(
    dat,
    yvar = yvar,
    tau = tau,
    h = h,
    method = method,
    df_spline = df_spline
  )$sigma_b2

  list(
    point = point,
    point_positive = max(0, point),
    boot = boot,
    boot_se = sd(boot),
    ci_percentile = unname(quantile(boot, c(0.025, 0.975))),
    p_two_sided_normal = 2 * pnorm(-abs(point / sd(boot))),
    p_upper_boot = mean(boot <= 0),
    B_ok = length(boot)
  )
}

data_path <- "../data/V1_V7_longformat_k1.csv"
dat_ftb <- prep_ftb(data_path)
B_BOOT <- as.integer(Sys.getenv("B_BOOT", "100"))
if (is.na(B_BOOT) || B_BOOT < 1) B_BOOT <- 100

cat("FTB rows:", nrow(dat_ftb), "\n")
cat("FTB subjects:", length(unique(dat_ftb$Participant_ID)), "\n\n")

quad <- estimate_sigma_b(dat_ftb, yvar = "SWS", tau = 0.5, method = "quadratic")
spline <- estimate_sigma_b(dat_ftb, yvar = "SWS", tau = 0.5, method = "spline")

cat("Smoke-test point estimates for SWS, tau = 0.5\n")
print(quad[c("sigma_b2", "sigma_b2_positive", "n_pairs", "n", "n_id")])
print(spline[c("sigma_b2", "sigma_b2_positive", "n_pairs", "n", "n_id")])

cat("\nCluster bootstrap, B =", B_BOOT, "\n")
quad_boot <- bootstrap_sigma_b(dat_ftb, yvar = "SWS", tau = 0.5, method = "quadratic", B = B_BOOT)
spline_boot <- bootstrap_sigma_b(dat_ftb, yvar = "SWS", tau = 0.5, method = "spline", B = B_BOOT)

print(list(quadratic = quad_boot[names(quad_boot) != "boot"]))
print(list(spline = spline_boot[names(spline_boot) != "boot"]))
