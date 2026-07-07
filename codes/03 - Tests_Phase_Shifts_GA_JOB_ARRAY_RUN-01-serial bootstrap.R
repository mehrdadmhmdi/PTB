#!/usr/bin/env Rscript

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SLURM ARRAY VERSION (one task = one biomarker-tau)
#
# Gestation-/conception-aligned phase-shift quantile spline analysis
# implementing Algorithm (A)-(E) using the prespecified working criterion described in the paper
#
# (A) Common-shift model:
#     Q_tau(Y_ij | t_ij, g_i) = alpha_tau + f_tau(t_ij + phi_tau g_i)
#     Estimate phi_hat by grid search over PHI_GRID minimizing check-loss.
#
# (B) Test H0_shift: phi_tau = 0
#     Use SUBJECT (cluster) bootstrap as written in the algorithm:
#       - resample subjects with replacement
#       - refit (A) to get T_shift* = Obj*(0) - min_phi Obj*(phi)
#       - p_shift = P(T_shift* >= T_shift_obs)
#     Also use the same bootstrap draws for percentile CI for phi.
#
# (C) Heterogeneous-shift extension:
#     Q_tau(Y_ij | t_ij, g_i, i) = alpha_tau + f_tau(t_ij + phi_i g_i)
#     where, for sPTB subjects, phi_i ~ N(phi_tau, sigma^2_{phi,tau}).
#
#     For each sigma in SIGMA_GRID (sigma is SD), fit by alternating updates:
#       - update spline coefficients given current {phi_i}
#       - update each subject-specific total shift phi_i
#       - update phi_tau = mean(phi_i) over sPTB subjects
#
#     Penalized objective:
#       Obj_het(sigma) = sum rho_tau(resid) + (1/(2 sigma^2)) sum (phi_i - phi_tau)^2
#
#     Select sigma_hat using a prespecified penalized-loss criterion:
#       - "pen_obj"       : crit = Obj_het(sigma)
#       - "ald_marginal"  : crit = Obj_het(sigma) + n1*log(sigma)
#
#     Report sigma^2_hat = sigma_hat^2.
#
# (D) Test H0_het: sigma^2_{phi,tau} = 0
#     Use SUBJECT-LEVEL wild (Rademacher) bootstrap as written:
#       - Fit common-shift at phi_hat to get fitted values and residuals
#       - Y* = Yhat + w_i * r_i (w_i constant within subject)
#       - Refit (A) and (C) on Y* and compute
#           T_het* = Obj_shift_hat* - Obj_het_hat*
#       - p_het = P(T_het* >= T_het_obs)
#
# (E) CI for sigma^2:
#     SUBJECT bootstrap: resample subjects, refit (A)+(C), record sigma_hat*^2,
#     take percentile CI.
#
# Output per task:
#   - phi_hat, phi CI, T_shift_obs, p_shift
#   - sigma_hat, sig2_hat, sig2 CI, T_het_obs, p_het
#   - diagnostics (optional but helpful)
#
# Notes on rq.fit.br():
#   rq.fit.br() takes a numeric design matrix x and response y.
#   If you want an intercept, include a column of 1's in x.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(quantreg)
  library(splines)
})

# -----------------------------
# User settings
# -----------------------------
vars2 <- c("SWS", "Midband", "Intercept", "AC", "Slope")
taus  <- c(0.2, 0.5, 0.8)

DF_SPLINE <- 3

# Phase-shift grids are in WEEKS.
PHI_GRID   <- seq(-8, 8, by = 0.5)   # common shift
PHI_I_GRID <- seq(-8, 8, by = 0.5)  # subject-specific total shifts phi_i

# sigma grid is SD in phi_i ~ N(phi_tau, sigma^2).
SIGMA_GRID <- c(0.25, 0.5, 1, 2, 4)

# Bootstrap counts
B_SHIFT <- 100   # subject bootstrap for (B): phi CI and p_shift (as in algorithm)
B_CI    <- 100   # subject bootstrap for (E): sigma^2 CI
B_HET   <- 100   # wild cluster bootstrap for (D): p_het

PTB_LEVEL <- 1
BASE_SEED <- 2026

MAX_ITER_RE <- 10
TOL_RE      <- 1e-3

# sigma selection criterion (prespecified):
#  - "pen_obj"      : crit = Obj_het(sigma)  [literal penalized loss]
#  - "ald_marginal" : crit = Obj_het(sigma) + n1*log(sigma)
#                    (ALD pseudo-marginal; prevents always selecting max sigma)
SIGMA_SELECT <- "ald_marginal"

# paths
DATA_PATH <- "/data/V1_V7_longformat_k1.csv"
OUT_DIR   <- "Test_results"
PART_DIR  <- file.path(OUT_DIR, "parts")
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(PART_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# SLURM array task selection
# -----------------------------
grid <- expand.grid(
  biomarker = vars2,
  tau = taus,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
if (is.na(task_id) || task_id < 1 || task_id > nrow(grid)) {
  stop(sprintf(
    "Invalid SLURM_ARRAY_TASK_ID=%s. Must be 1..%d",
    Sys.getenv("SLURM_ARRAY_TASK_ID"), nrow(grid)
  ))
}

yvar <- grid$biomarker[task_id]
tau  <- grid$tau[task_id]

tau_tag   <- gsub("\\.", "p", format(tau, trim = TRUE))
part_file <- file.path(
  PART_DIR,
  sprintf("qr_phase_shift_%s_tau%s_task%02d.csv", yvar, tau_tag, task_id)
)

if (file.exists(part_file)) {
  message("Part file exists, skipping: ", part_file)
  quit(save = "no", status = 0)
}

message("TASK ", task_id, "/", nrow(grid), "  biomarker=", yvar, "  tau=", tau)

# -----------------------------
# Data load + prep
# -----------------------------
resolve_ptb_level <- function(dat, outcome_col, ptb_level) {
  levs <- levels(dat[[outcome_col]])
  if (is.numeric(ptb_level) && length(ptb_level) == 1L) {
    if (ptb_level >= 1 && ptb_level <= length(levs)) return(levs[ptb_level])
  }
  if (is.character(ptb_level) && ptb_level %in% levs) return(ptb_level)
  stop("ptb_level must be a valid Outcome level name or a valid factor level index.")
}

raw_dat <- read_csv(DATA_PATH, show_col_types = FALSE)

# IMPORTANT: GA scale guard.
# If GA looks like proportion (<= ~2), convert to weeks by *39.
# If GA already looks like weeks, do NOT multiply.
ga_raw <- suppressWarnings(as.numeric(raw_dat$GA))
ga_max <- max(ga_raw, na.rm = TRUE)

ga_mult <- if (is.finite(ga_max) && ga_max <= 2) 39 else 1
message("GA scaling: max(raw GA)=", signif(ga_max, 4), "  -> multiplying by ", ga_mult)

df0 <- raw_dat %>%
  mutate(
    GA = as.numeric(GA) * ga_mult,
    Outcome = as.factor(Outcome),
    Participant_ID = as.character(Participant_ID)
  )

ptb_lab <- resolve_ptb_level(df0, "Outcome", PTB_LEVEL)

df <- df0 %>%
  mutate(G = as.integer(Outcome == ptb_lab)) %>%
  filter(
    !is.na(.data[[yvar]]),
    !is.na(GA),
    !is.na(G),
    !is.na(Participant_ID)
  ) %>%
  arrange(Participant_ID, GA)

if (nrow(df) == 0) stop("No usable rows after filtering.")
if (length(unique(df$Participant_ID)) < 4) stop("Too few subjects after filtering.")

# Sanity: ensure G constant within subject
if (any(tapply(df$G, df$Participant_ID, function(z) any(z != z[1])))) {
  stop("G is not constant within Participant_ID. Fix data first.")
}

# -----------------------------
# Core helpers
# -----------------------------
qr_checkloss <- function(u, tau) {
  # sum rho_tau(u)
  sum(u * (tau - (u < 0)))
}

safe_quantile <- function(x, probs) {
  if (all(is.na(x))) return(rep(NA_real_, length(probs)))
  as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

build_spline_spec <- function(x_base, df_spline, lower_extra = 0, upper_extra = 0) {
  # df_spline here controls the spline basis dimension (including boundary behavior).
  n_int <- max(df_spline - 1L, 0L)
  
  knots <- NULL
  if (n_int > 0) {
    probs <- seq(0, 1, length.out = n_int + 2L)[-c(1, n_int + 2L)]
    knots <- as.numeric(quantile(x_base, probs = probs, na.rm = TRUE, type = 7))
    knots <- unique(knots)
    if (length(knots) == 0) knots <- NULL
  }
  
  bknots <- c(
    min(x_base, na.rm = TRUE) + lower_extra,
    max(x_base, na.rm = TRUE) + upper_extra
  )
  
  list(knots = knots, Boundary.knots = bknots)
}

make_X <- function(x_shift, spline_spec) {
  xb <- ns(
    x_shift,
    knots = spline_spec$knots,
    Boundary.knots = spline_spec$Boundary.knots,
    intercept = FALSE
  )
  X <- cbind(`(Intercept)` = 1, xb)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  X
}

fit_qr <- function(X, y, tau) {
  fit <- tryCatch(
    rq.fit.br(x = X, y = y, tau = tau),
    error = function(e) NULL
  )
  if (is.null(fit) || is.null(fit$coefficients) || anyNA(fit$coefficients)) return(NULL)
  
  beta <- as.numeric(fit$coefficients)
  resid <- y - drop(X %*% beta)
  list(beta = beta, resid = resid, obj = qr_checkloss(resid, tau))
}

fit_fixed_phi <- function(dat, yvar, tau, spline_spec, phi,
                          t_col = "GA", g_col = "G") {
  y <- dat[[yvar]]
  t <- dat[[t_col]]
  g <- dat[[g_col]]
  
  x_shift <- t + phi * g
  X <- make_X(x_shift, spline_spec)
  fit <- fit_qr(X, y, tau)
  
  if (is.null(fit)) {
    return(list(phi = phi, obj = Inf, beta = NULL,
                fitted = rep(NA_real_, length(y)),
                resid  = rep(NA_real_, length(y))))
  }
  
  fitted <- y - fit$resid
  list(phi = phi, obj = fit$obj, beta = fit$beta, fitted = fitted, resid = fit$resid)
}

# Step (A): profile over PHI_GRID, choose phi_hat = argmin Obj(phi)
fit_common_shift <- function(dat, yvar, tau, spline_spec, phi_grid,
                             t_col = "GA", g_col = "G") {
  objs <- rep(Inf, length(phi_grid))
  for (k in seq_along(phi_grid)) {
    objs[k] <- fit_fixed_phi(dat, yvar, tau, spline_spec, phi_grid[k], t_col, g_col)$obj
  }
  if (!any(is.finite(objs))) {
    return(list(phi_hat = NA_real_, obj_min = NA_real_, obj0 = NA_real_,
                T_shift = NA_real_, fitted_hat = NULL, resid_hat = NULL,
                profile = data.frame(phi = phi_grid, obj = objs)))
  }
  
  k_min <- which.min(objs)
  phi_hat <- phi_grid[k_min]
  obj_min <- objs[k_min]
  
  fit_hat <- fit_fixed_phi(dat, yvar, tau, spline_spec, phi_hat, t_col, g_col)
  fit_0   <- fit_fixed_phi(dat, yvar, tau, spline_spec, 0,      t_col, g_col)
  
  list(
    phi_hat = phi_hat,
    obj_min = obj_min,
    obj0 = fit_0$obj,
    T_shift = fit_0$obj - obj_min,
    fitted_hat = fit_hat$fitted,
    resid_hat  = fit_hat$resid,
    profile = data.frame(phi = phi_grid, obj = objs)
  )
}

# Subject bootstrap (cluster resampling)
resample_subjects <- function(dat, id_col = "Participant_ID") {
  ids <- unique(dat[[id_col]])
  n_ids <- length(ids)
  samp_ids <- sample(ids, size = n_ids, replace = TRUE)
  
  out <- vector("list", length(samp_ids))
  for (b in seq_along(samp_ids)) {
    tmp <- dat[dat[[id_col]] == samp_ids[b], , drop = FALSE]
    tmp[[id_col]] <- paste0("boot_", b)  # new IDs, preserve clustering
    out[[b]] <- tmp
  }
  bind_rows(out)
}

# -----------------------------
# Step (C): Heterogeneous fit for a fixed sigma
# -----------------------------
fit_het_given_sigma <- function(dat, yvar, tau, spline_spec,
                                sigma,
                                phi_i_grid,
                                id_col = "Participant_ID",
                                t_col = "GA",
                                g_col = "G",
                                init_phi,
                                max_iter = 10,
                                tol = 1e-3) {
  y  <- dat[[yvar]]
  t  <- dat[[t_col]]
  g  <- dat[[g_col]]
  id <- dat[[id_col]]
  
  ptb_ids <- unique(id[g == 1L])
  n1 <- length(ptb_ids)
  if (n1 == 0L) return(NULL)
  
  # Initialize subject-specific total shifts at the common-shift estimate
  phi <- init_phi
  phi_i_by_id <- setNames(rep(init_phi, n1), ptb_ids)
  
  phi_vec_from_map <- function(id_vec, g_vec, phi_map) {
    out <- rep(0, length(id_vec))
    idx <- which(g_vec == 1L)
    if (length(idx) > 0) {
      out[idx] <- unname(phi_map[id_vec[idx]])
      out[is.na(out)] <- 0
    }
    out
  }
  
  full_fit <- function(phi_val, phi_map) {
    phi_i_vec <- phi_vec_from_map(id, g, phi_map)
    x_shift <- t + phi_i_vec * g
    X <- make_X(x_shift, spline_spec)
    fit <- fit_qr(X, y, tau)
    if (is.null(fit)) return(NULL)
    
    pen <- sum((phi_map - phi_val)^2) / (2 * sigma^2)
    
    list(
      beta = fit$beta,
      check_obj = fit$obj,
      pen_obj = fit$obj + pen,
      phi_i_vec = phi_i_vec
    )
  }
  
  for (it in seq_len(max_iter)) {
    phi_old   <- phi
    phi_i_old <- phi_i_by_id
    
    # (1) Fit spline coefficients given current {phi_i}
    cur <- full_fit(phi, phi_i_by_id)
    if (is.null(cur)) return(NULL)
    beta_cur <- cur$beta
    
    # (2) Update each subject-specific total shift phi_i
    for (sid in ptb_ids) {
      idx_i <- which(id == sid & g == 1L)
      if (length(idx_i) == 0L) next
      
      t_i <- t[idx_i]
      y_i <- y[idx_i]
      
      loss <- rep(Inf, length(phi_i_grid))
      for (kk in seq_along(phi_i_grid)) {
        phi_i_cand <- phi_i_grid[kk]
        X_i <- make_X(t_i + phi_i_cand, spline_spec)
        u_i <- y_i - drop(X_i %*% beta_cur)
        loss[kk] <- qr_checkloss(u_i, tau) + (phi_i_cand - phi)^2 / (2 * sigma^2)
      }
      phi_i_by_id[sid] <- phi_i_grid[which.min(loss)]
    }
    
    # (3) Update population-average shift
    phi <- mean(phi_i_by_id)
    
    dphi   <- abs(phi - phi_old)
    dphi_i <- max(abs(phi_i_by_id - phi_i_old))
    if (max(dphi, dphi_i) < tol) break
  }
  
  fin <- full_fit(phi, phi_i_by_id)
  if (is.null(fin)) return(NULL)
  
  list(
    phi_hat = phi,
    phi_i_hat = phi_i_by_id,
    sigma_hat = sigma,
    sigma2_hat = sigma^2,
    n1 = n1,
    check_obj = fin$check_obj,
    pen_obj = fin$pen_obj
  )
}

# Select sigma_hat by prespecified criterion; also returns Obj_het_hat and the fit
fit_het_select_sigma <- function(dat, yvar, tau, spline_spec,
                                 phi_i_grid, sigma_grid,
                                 id_col = "Participant_ID",
                                 t_col = "GA",
                                 g_col = "G",
                                 init_phi,
                                 max_iter = 10,
                                 tol = 1e-3,
                                 sigma_select = "ald_marginal") {
  g_first <- tapply(dat[[g_col]], dat[[id_col]], function(z) z[1])
  n1 <- sum(g_first == 1L, na.rm = TRUE)
  if (n1 <= 1L) return(NULL)
  
  best_fit <- NULL
  best_crit <- Inf
  
  for (sg in sigma_grid) {
    fit <- fit_het_given_sigma(
      dat = dat, yvar = yvar, tau = tau, spline_spec = spline_spec,
      sigma = sg, phi_i_grid = phi_i_grid,
      id_col = id_col, t_col = t_col, g_col = g_col,
      init_phi = init_phi,
      max_iter = max_iter, tol = tol
    )
    if (is.null(fit) || !is.finite(fit$pen_obj)) next
    
    crit <- switch(
      sigma_select,
      pen_obj = fit$pen_obj,
      ald_marginal = fit$pen_obj + n1 * log(sg),
      stop("Unknown sigma_select. Use 'pen_obj' or 'ald_marginal'.")
    )
    
    if (crit < best_crit) {
      best_crit <- crit
      best_fit <- fit
      best_fit$crit <- crit
    }
  }
  
  best_fit
}

# -----------------------------
# Full analysis (A)-(E) for one biomarker/tau
# -----------------------------
analyze_one_task <- function(dat, yvar, tau,
                             id_col = "Participant_ID",
                             t_col = "GA",
                             g_col = "G",
                             df_spline = 3,
                             phi_grid = PHI_GRID,
                             phi_i_grid = PHI_I_GRID,
                             sigma_grid = SIGMA_GRID,
                             B_shift = B_SHIFT,
                             B_ci    = B_CI,
                             B_het   = B_HET,
                             max_iter_re = MAX_ITER_RE,
                             tol_re = TOL_RE,
                             sigma_select = SIGMA_SELECT,
                             seed = BASE_SEED) {
  set.seed(seed)
  
  # Fix spline basis once (Algorithm uses fixed B(.))
  lower_extra <- min(c(0, phi_grid, phi_i_grid), na.rm = TRUE)
  upper_extra <- max(c(0, phi_grid, phi_i_grid), na.rm = TRUE)
  
  spline_spec <- build_spline_spec(
    x_base = dat[[t_col]],
    df_spline = df_spline,
    lower_extra = lower_extra,
    upper_extra = upper_extra
  )
  
  n_obs <- nrow(dat)
  n_id  <- length(unique(dat[[id_col]]))
  n_ptb <- length(unique(dat[[id_col]][dat[[g_col]] == 1L]))
  n_ftb <- length(unique(dat[[id_col]][dat[[g_col]] == 0L]))
  
  # ============================================================
  # (A) Common-shift model
  # ============================================================
  common_fit <- fit_common_shift(dat, yvar, tau, spline_spec, phi_grid, t_col, g_col)
  if (!is.finite(common_fit$T_shift)) stop("Common-shift model failed.")
  
  phi_hat <- common_fit$phi_hat
  Obj_shift_hat <- common_fit$obj_min      # \widehat{Obj}^{shift}
  Obj0 <- common_fit$obj0                  # Obj(0)
  T_shift_obs <- common_fit$T_shift        # Obj(0) - min_phi Obj(phi)
  
  # ============================================================
  # (B) Shift bootstrap for CI and p-value (subject bootstrap)
  # ============================================================
  set.seed(seed + 100)
  phi_boot <- rep(NA_real_, B_shift)
  T_shift_boot <- rep(NA_real_, B_shift)
  
  for (b in seq_len(B_shift)) {
    dat_b <- resample_subjects(dat, id_col = id_col)
    fit_b <- fit_common_shift(dat_b, yvar, tau, spline_spec, phi_grid, t_col, g_col)
    phi_boot[b] <- fit_b$phi_hat
    T_shift_boot[b] <- fit_b$T_shift
  }
  
  phi_ci <- safe_quantile(phi_boot, c(0.025, 0.975))
  
  # Algorithm says "compute bootstrap p-value using {T_shift*}".
  # We implement p = P(T* >= T_obs) using the bootstrap draws.
  p_shift <- mean(T_shift_boot >= T_shift_obs, na.rm = TRUE)
  
  # ============================================================
  # (C) Heterogeneous model: fit and select sigma_hat
  # ============================================================
  sigma_hat <- NA_real_
  sig2_hat  <- NA_real_
  Obj_het_hat <- NA_real_
  Crit_het_hat <- NA_real_
  T_het_obs <- NA_real_
  phi_i_var_hat <- NA_real_
  n_nonzero_dev_phi_i <- NA_integer_
  max_abs_dev_phi_i <- NA_real_
  
  het_fit <- NULL
  if (n_ptb > 1L) {
    het_fit <- fit_het_select_sigma(
      dat = dat, yvar = yvar, tau = tau, spline_spec = spline_spec,
      phi_i_grid = phi_i_grid, sigma_grid = sigma_grid,
      id_col = id_col, t_col = t_col, g_col = g_col,
      init_phi = phi_hat,
      max_iter = max_iter_re, tol = tol_re,
      sigma_select = sigma_select
    )
    
    if (!is.null(het_fit) && is.finite(het_fit$pen_obj)) {
      sigma_hat  <- het_fit$sigma_hat
      sig2_hat   <- het_fit$sigma2_hat
      Obj_het_hat <- het_fit$pen_obj
      Crit_het_hat <- het_fit$crit
      
      # Diagnostics on subject-specific shifts phi_i
      phi_i_vec <- as.numeric(het_fit$phi_i_hat)
      dev_vec <- phi_i_vec - het_fit$phi_hat
      phi_i_var_hat <- if (length(phi_i_vec) > 1) var(phi_i_vec) else 0
      n_nonzero_dev_phi_i <- sum(abs(dev_vec) > 1e-8)
      max_abs_dev_phi_i <- if (length(dev_vec) > 0) max(abs(dev_vec)) else 0
      
      # ============================================================
      # (D) Heterogeneity test statistic (algorithm definition)
      # ============================================================
      # T_het = Obj_shift_hat - Obj_het_hat  (penalized het objective)
      T_het_obs <- Obj_shift_hat - Obj_het_hat
    }
  }
  
  # ============================================================
  # (D) Wild cluster bootstrap under H0_het
  # ============================================================
  p_het <- NA_real_
  T_het_boot <- rep(NA_real_, B_het)
  
  if (is.finite(T_het_obs) && n_ptb > 1L) {
    # Fit common shift at phi_hat, get fitted values and residuals
    fit_shift_hat <- fit_fixed_phi(dat, yvar, tau, spline_spec, phi_hat, t_col, g_col)
    if (!all(is.finite(fit_shift_hat$fitted)) || !all(is.finite(fit_shift_hat$resid))) {
      warning("Common-shift fitted/resid contain NA; cannot do heterogeneity bootstrap reliably.")
    } else {
      set.seed(seed + 200)
      id_vals <- unique(dat[[id_col]])
      id_row  <- dat[[id_col]]
      
      for (b in seq_len(B_het)) {
        w_id <- setNames(sample(c(-1, 1), size = length(id_vals), replace = TRUE), id_vals)
        
        dat_star <- dat
        dat_star[[yvar]] <- fit_shift_hat$fitted + unname(w_id[id_row]) * fit_shift_hat$resid
        
        # Refit (A)
        common_star <- fit_common_shift(dat_star, yvar, tau, spline_spec, phi_grid, t_col, g_col)
        if (!is.finite(common_star$obj_min)) next
        
        # Refit (C) and select sigma
        het_star <- fit_het_select_sigma(
          dat = dat_star, yvar = yvar, tau = tau, spline_spec = spline_spec,
          phi_i_grid = phi_i_grid, sigma_grid = sigma_grid,
          id_col = id_col, t_col = t_col, g_col = g_col,
          init_phi = common_star$phi_hat,
          max_iter = max_iter_re, tol = tol_re,
          sigma_select = sigma_select
        )
        if (is.null(het_star) || !is.finite(het_star$pen_obj)) next
        
        # Compute T_het* with the SAME definition:
        # Obj_shift_hat* - Obj_het_hat* (penalized het objective)
        T_het_boot[b] <- common_star$obj_min - het_star$pen_obj
      }
      
      p_het <- mean(T_het_boot >= T_het_obs, na.rm = TRUE)
    }
  }
  
  # ============================================================
  # (E) CI for sigma^2 via subject bootstrap
  # ============================================================
  sig2_ci <- c(NA_real_, NA_real_)
  if (n_ptb > 1L) {
    set.seed(seed + 300)
    sig2_boot <- rep(NA_real_, B_ci)
    
    for (b in seq_len(B_ci)) {
      dat_b <- resample_subjects(dat, id_col = id_col)
      
      # Refit (A)
      common_b <- fit_common_shift(dat_b, yvar, tau, spline_spec, phi_grid, t_col, g_col)
      if (!is.finite(common_b$phi_hat)) next
      
      # Refit (C) and select sigma
      het_b <- fit_het_select_sigma(
        dat = dat_b, yvar = yvar, tau = tau, spline_spec = spline_spec,
        phi_i_grid = phi_i_grid, sigma_grid = sigma_grid,
        id_col = id_col, t_col = t_col, g_col = g_col,
        init_phi = common_b$phi_hat,
        max_iter = max_iter_re, tol = tol_re,
        sigma_select = sigma_select
      )
      sig2_boot[b] <- if (is.null(het_b)) NA_real_ else het_b$sigma2_hat
    }
    
    sig2_ci <- safe_quantile(sig2_boot, c(0.025, 0.975))
  }
  
  # Return results
  list(
    biomarker = yvar,
    tau = tau,
    n_obs = n_obs,
    n_id = n_id,
    n_ptb = n_ptb,
    n_ftb = n_ftb,
    
    # (A)-(B) common-shift
    phi_hat = phi_hat,
    phi_ci_low = phi_ci[1],
    phi_ci_high = phi_ci[2],
    Obj_shift_hat = Obj_shift_hat,
    Obj0 = Obj0,
    T_shift_obs = T_shift_obs,
    p_shift = p_shift,
    
    # (C)-(E) heterogeneity
    sigma_hat = sigma_hat,
    sig2_hat = sig2_hat,
    sig2_ci_low = sig2_ci[1],
    sig2_ci_high = sig2_ci[2],
    Obj_het_hat = Obj_het_hat,
    Crit_het_hat = Crit_het_hat,
    T_het_obs = T_het_obs,
    p_het = p_het,
    
    # diagnostics
    phi_i_var_hat = phi_i_var_hat,
    n_nonzero_dev_phi_i = n_nonzero_dev_phi_i,
    max_abs_dev_phi_i = max_abs_dev_phi_i,
    sigma_select = sigma_select
  )
}

# -----------------------------
# Run one task
# -----------------------------
res <- tryCatch(
  analyze_one_task(
    dat = df,
    yvar = yvar,
    tau = tau,
    df_spline = DF_SPLINE,
    phi_grid = PHI_GRID,
    phi_i_grid = PHI_I_GRID,
    sigma_grid = SIGMA_GRID,
    B_shift = B_SHIFT,
    B_ci = B_CI,
    B_het = B_HET,
    max_iter_re = MAX_ITER_RE,
    tol_re = TOL_RE,
    sigma_select = SIGMA_SELECT,
    seed = BASE_SEED + task_id
  ),
  error = function(e) {
    list(biomarker = yvar, tau = tau, error = conditionMessage(e))
  }
)

out1 <- as.data.frame(res, stringsAsFactors = FALSE)
write.csv(out1, part_file, row.names = FALSE)
message("Wrote: ", part_file)
