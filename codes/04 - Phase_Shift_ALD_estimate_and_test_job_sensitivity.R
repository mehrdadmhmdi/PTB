#!/usr/bin/env Rscript

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SLURM ARRAY VERSION (one task = one biomarker-tau)
#
# Conception-aligned phase-shift analysis with sensitivity to the random-effects
# distribution in a random-shift ALD mixed model:
#
#   Y_ij = alpha_tau + B(t_ij + phi_i g_i)^T beta_tau + eps_ij,tau
#   eps_ij,tau ~ ALD_tau(0, s_tau)
#
#   phi_i = phi_tau + sigma_phi,tau * U_i
#   U_i ~ G, where G is varied over several parametric specifications
#         (Normal, Student-t with fixed df, Laplace)
#
# Key implementation details:
#   1) Keep the OUTER SLURM array over biomarker–tau combinations.
#   2) Precompute all fixed spline objects and common-model QR starts.
#      Exact full precomputation of all shifted spline matrices is impossible
#      because phi and sigma are continuous, but all fixed objects are cached.
#   3) Parallelize ONLY the bootstrap replicates inside each task.
#
# Output:
#   One CSV per biomarker-tau task with one row per random-effects distribution.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

# Avoid BLAS/OpenMP oversubscription when bootstrap replicates are parallelized.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(quantreg)
  library(splines)
  library(parallel)
})

# -----------------------------
# User settings
# -----------------------------
vars2 <- c("SWS", "Midband", "Intercept", "AC", "Slope")
taus  <- c(0.2, 0.5, 0.8)

DF_SPLINE <- 3

# Common-shift initialization / plausible bounds
PHI_GRID <- seq(-8, 8, by = 0.5)
PHI_LOWER <- min(PHI_GRID)
PHI_UPPER <- max(PHI_GRID)

# Random-effect sensitivity specifications
# family in {"normal", "student_t", "laplace"}
DIST_SPECS <- list(
  list(name = "normal",         family = "normal",    param = NA_real_),
  list(name = "student_t_df4",  family = "student_t", param = 4),
  list(name = "student_t_df8",  family = "student_t", param = 8),
  list(name = "laplace",        family = "laplace",   param = NA_real_)
)

NQ_DIST <- 31  # quantile quadrature nodes for the random effect

# Starts / bounds
SIGMA_INIT_GRID <- c(0.15, 0.30, 0.60, 1.00)
SIGMA_LOWER <- 1e-4
SIGMA_UPPER <- 8
SCALE_LOWER <- 1e-5
SCALE_UPPER <- 1e3

PTB_LEVEL <- 1
BASE_SEED <- 2026

# Bootstrap counts
B_SHIFT_COMMON_TEST <- 100   # common-model test of phi=0
B_SHIFT_COMMON_CI   <- 100   # common-model phi CI

B_SHIFT_HET_TEST <- 100      # hetero-model test of phi=0 (spec-specific)
B_SHIFT_HET_CI   <- 100      # hetero-model phi CI (spec-specific)
B_HET_TEST       <- 100      # heterogeneity test sigma^2=0 (spec-specific)
B_SIG2_CI        <- 100      # sigma^2 CI (spec-specific)

# Optimization controls
OPT_MAXIT_COMMON <- 300
OPT_MAXIT_HET    <- 300

# Number of worker processes inside each SLURM task
NCORES_TASK <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))

# paths
DATA_PATH <- "data/V1_V7_longformat_k1.csv"
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
  sprintf("ald_sensitivity_%s_tau%s_task%02d.csv", yvar, tau_tag, task_id)
)

if (file.exists(part_file)) {
  message("Part file exists, skipping: ", part_file)
  quit(save = "no", status = 0)
}

message(
  "TASK ", task_id, "/", nrow(grid),
  "  biomarker=", yvar,
  "  tau=", tau,
  "  cpus-per-task=", NCORES_TASK
)

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

if (any(tapply(df$G, df$Participant_ID, function(z) any(z != z[1])))) {
  stop("G is not constant within Participant_ID. Fix data first.")
}

# -----------------------------
# Core helpers
# -----------------------------
rho_tau_vec <- function(u, tau) {
  u * (tau - (u < 0))
}

rho_tau_sum <- function(u, tau) {
  sum(rho_tau_vec(u, tau))
}

safe_quantile <- function(x, probs) {
  if (all(is.na(x))) return(rep(NA_real_, length(probs)))
  as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

log_sum_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

ald_logpdf <- function(y, mu, tau, s) {
  log(tau * (1 - tau)) - log(s) - rho_tau_vec(y - mu, tau) / s
}

rALD <- function(n, tau, s) {
  u <- runif(n)
  x <- numeric(n)
  left <- (u < tau)
  x[left]  <- s / (1 - tau) * log(u[left] / tau)
  x[!left] <- -s / tau * log((1 - u[!left]) / (1 - tau))
  x
}

safe_name_num <- function(x, digits = 3) {
  if (is.na(x)) return(NA_character_)
  formatC(x, digits = digits, format = "fg", flag = "#")
}

# -----------------------------
# Random-effect distributions
# -----------------------------
std_t_scale <- function(nu) sqrt((nu - 2) / nu)

qstdt <- function(p, nu) std_t_scale(nu) * qt(p, df = nu)
rstdt <- function(n, nu) std_t_scale(nu) * rt(n, df = nu)

# Standardized Laplace: mean 0, variance 1
laplace_b <- 1 / sqrt(2)

qstdlaplace <- function(p) {
  b <- laplace_b
  out <- numeric(length(p))
  left <- (p < 0.5)
  out[left]  <- b * log(2 * p[left])
  out[!left] <- -b * log(2 * (1 - p[!left]))
  out
}

rstdlaplace <- function(n) {
  qstdlaplace(runif(n))
}

build_quad_rule <- function(dist_spec, nq = 31) {
  p <- (seq_len(nq) - 0.5) / nq
  
  nodes <- switch(
    dist_spec$family,
    normal    = qnorm(p),
    student_t = qstdt(p, nu = dist_spec$param),
    laplace   = qstdlaplace(p),
    stop("Unknown random-effect family.")
  )
  
  list(
    nodes = nodes,
    weights = rep(1 / nq, nq),
    nq = nq,
    dist_name = dist_spec$name,
    family = dist_spec$family,
    param = dist_spec$param
  )
}

r_random_effect_std <- function(n, dist_spec) {
  switch(
    dist_spec$family,
    normal    = rnorm(n),
    student_t = rstdt(n, nu = dist_spec$param),
    laplace   = rstdlaplace(n),
    stop("Unknown random-effect family.")
  )
}

# -----------------------------
# Spline helpers
# -----------------------------
build_spline_spec <- function(x_base, df_spline, lower_extra = 0, upper_extra = 0) {
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

make_B <- function(x_shift, spline_spec) {
  B <- ns(
    x_shift,
    knots = spline_spec$knots,
    Boundary.knots = spline_spec$Boundary.knots,
    intercept = FALSE
  )
  B <- as.matrix(B)
  storage.mode(B) <- "double"
  B
}

make_X_qr <- function(x_shift, spline_spec) {
  B <- make_B(x_shift, spline_spec)
  X <- cbind(`(Intercept)` = 1, B)
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
  
  beta  <- as.numeric(fit$coefficients)
  resid <- y - drop(X %*% beta)
  list(beta = beta, resid = resid, obj = rho_tau_sum(resid, tau))
}

parallel_lapply_safe <- function(X, FUN, mc.cores = 1L) {
  mc.cores <- max(1L, as.integer(mc.cores))
  if (.Platform$OS.type == "windows" || mc.cores <= 1L || length(X) <= 1L) {
    lapply(X, FUN)
  } else {
    mclapply(
      X,
      FUN,
      mc.cores = min(mc.cores, length(X)),
      mc.preschedule = FALSE,
      mc.set.seed = FALSE
    )
  }
}

# -----------------------------
# Precomputed base context
# -----------------------------
build_base_context <- function(dat, yvar, spline_spec, phi_grid,
                               id_col = "Participant_ID",
                               t_col = "GA",
                               g_col = "G") {
  t  <- dat[[t_col]]
  g  <- dat[[g_col]]
  id <- dat[[id_col]]
  
  subj_ids <- unique(id)
  subj_index <- setNames(lapply(subj_ids, function(s) which(id == s)), subj_ids)
  subj_g <- sapply(subj_ids, function(s) g[subj_index[[s]][1]])
  
  ptb_ids <- subj_ids[subj_g == 1L]
  ftb_ids <- subj_ids[subj_g == 0L]
  
  ptb_idx <- which(g == 1L)
  ftb_idx <- which(g == 0L)
  
  ptb_subject_factor <- factor(id[ptb_idx], levels = ptb_ids)
  
  B0_full <- make_B(t, spline_spec)
  
  idx_phi_zero <- match(0, phi_grid)
  X_common_qr <- vector("list", length(phi_grid))
  for (k in seq_along(phi_grid)) {
    phi <- phi_grid[k]
    X_common_qr[[k]] <- make_X_qr(t + phi * g, spline_spec)
  }
  
  list(
    yvar = yvar,
    t = t,
    g = g,
    id = id,
    spline_spec = spline_spec,
    phi_grid = phi_grid,
    idx_phi_zero = idx_phi_zero,
    B0_full = B0_full,
    X_common_qr = X_common_qr,
    subj_ids = subj_ids,
    subj_index = subj_index,
    subj_g = subj_g,
    ptb_ids = ptb_ids,
    ftb_ids = ftb_ids,
    ptb_idx = ptb_idx,
    ftb_idx = ftb_idx,
    ptb_subject_factor = ptb_subject_factor,
    n_basis = ncol(B0_full)
  )
}

# -----------------------------
# Common-model QR initialization
# -----------------------------
fit_common_qr_init <- function(base_ctx, y, tau, phi_fixed = NULL) {
  if (!is.null(phi_fixed)) {
    X <- if (abs(phi_fixed) < .Machine$double.eps^0.5 && !is.na(base_ctx$idx_phi_zero)) {
      base_ctx$X_common_qr[[base_ctx$idx_phi_zero]]
    } else {
      make_X_qr(base_ctx$t + phi_fixed * base_ctx$g, base_ctx$spline_spec)
    }
    fit <- fit_qr(X, y, tau)
    if (is.null(fit)) {
      alpha0 <- as.numeric(quantile(y, probs = tau, na.rm = TRUE, names = FALSE))
      beta0 <- rep(0, base_ctx$n_basis)
      s0 <- max(mean(rho_tau_vec(y - alpha0, tau)), 1e-3)
      return(list(alpha = alpha0, beta = beta0, phi = phi_fixed, s = s0))
    }
    alpha0 <- fit$beta[1]
    beta0  <- fit$beta[-1]
    s0 <- max(mean(rho_tau_vec(fit$resid, tau)), 1e-3)
    return(list(alpha = alpha0, beta = beta0, phi = phi_fixed, s = s0))
  }
  
  objs <- rep(Inf, length(base_ctx$phi_grid))
  fits <- vector("list", length(base_ctx$phi_grid))
  for (k in seq_along(base_ctx$phi_grid)) {
    fits[[k]] <- fit_qr(base_ctx$X_common_qr[[k]], y, tau)
    if (!is.null(fits[[k]])) objs[k] <- fits[[k]]$obj
  }
  
  if (!any(is.finite(objs))) {
    alpha0 <- as.numeric(quantile(y, probs = tau, na.rm = TRUE, names = FALSE))
    beta0 <- rep(0, base_ctx$n_basis)
    s0 <- max(mean(rho_tau_vec(y - alpha0, tau)), 1e-3)
    return(list(alpha = alpha0, beta = beta0, phi = 0, s = s0))
  }
  
  kmin <- which.min(objs)
  fit <- fits[[kmin]]
  alpha0 <- fit$beta[1]
  beta0  <- fit$beta[-1]
  phi0   <- base_ctx$phi_grid[kmin]
  s0 <- max(mean(rho_tau_vec(fit$resid, tau)), 1e-3)
  
  list(alpha = alpha0, beta = beta0, phi = phi0, s = s0)
}

# -----------------------------
# Common-shift ALD MLE
# -----------------------------
common_negloglik <- function(par, base_ctx, y, tau, phi_fixed = NULL) {
  k <- base_ctx$n_basis
  
  if (is.null(phi_fixed)) {
    alpha <- par[1]
    beta  <- par[2:(1 + k)]
    phi   <- par[2 + k]
    s     <- exp(par[3 + k])
  } else {
    alpha <- par[1]
    beta  <- par[2:(1 + k)]
    phi   <- phi_fixed
    s     <- exp(par[2 + k])
  }
  
  B <- make_B(base_ctx$t + phi * base_ctx$g, base_ctx$spline_spec)
  mu <- alpha + drop(B %*% beta)
  
  nll <- -sum(ald_logpdf(y, mu, tau, s))
  if (!is.finite(nll)) nll <- 1e100
  nll
}

fit_common_mle <- function(base_ctx, y, tau, phi_fixed = NULL,
                           start_fit = NULL,
                           maxit = OPT_MAXIT_COMMON) {
  k <- base_ctx$n_basis
  
  if (is.null(start_fit)) {
    init <- fit_common_qr_init(base_ctx, y, tau, phi_fixed = phi_fixed)
  } else {
    init <- list(
      alpha = start_fit$alpha,
      beta  = start_fit$beta,
      phi   = if (is.null(phi_fixed)) start_fit$phi else phi_fixed,
      s     = start_fit$s
    )
  }
  
  if (is.null(phi_fixed)) {
    par0 <- c(init$alpha, init$beta, init$phi, log(init$s))
    lower <- c(rep(-Inf, 1 + k), PHI_LOWER, log(SCALE_LOWER))
    upper <- c(rep( Inf, 1 + k), PHI_UPPER, log(SCALE_UPPER))
  } else {
    par0 <- c(init$alpha, init$beta, log(init$s))
    lower <- c(rep(-Inf, 1 + k), log(SCALE_LOWER))
    upper <- c(rep( Inf, 1 + k), log(SCALE_UPPER))
  }
  
  opt <- tryCatch(
    optim(
      par = par0,
      fn = common_negloglik,
      base_ctx = base_ctx,
      y = y,
      tau = tau,
      phi_fixed = phi_fixed,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = maxit)
    ),
    error = function(e) NULL
  )
  if (is.null(opt) || !is.finite(opt$value)) return(NULL)
  
  if (is.null(phi_fixed)) {
    alpha <- opt$par[1]
    beta  <- opt$par[2:(1 + k)]
    phi   <- opt$par[2 + k]
    s     <- exp(opt$par[3 + k])
  } else {
    alpha <- opt$par[1]
    beta  <- opt$par[2:(1 + k)]
    phi   <- phi_fixed
    s     <- exp(opt$par[2 + k])
  }
  
  B <- make_B(base_ctx$t + phi * base_ctx$g, base_ctx$spline_spec)
  mu <- alpha + drop(B %*% beta)
  loglik <- -common_negloglik(opt$par, base_ctx, y, tau, phi_fixed = phi_fixed)
  
  list(
    alpha = alpha,
    beta = beta,
    phi = phi,
    s = s,
    mu = mu,
    resid = y - mu,
    loglik = loglik,
    conv = opt$convergence
  )
}

# -----------------------------
# Heterogeneous random-shift ALD MLE
# -----------------------------
hetero_negloglik <- function(par, base_ctx, y, tau, quad_rule, phi_fixed = NULL) {
  k <- base_ctx$n_basis
  
  if (is.null(phi_fixed)) {
    alpha <- par[1]
    beta  <- par[2:(1 + k)]
    phi   <- par[2 + k]
    sigma <- exp(par[3 + k])
    s     <- exp(par[4 + k])
  } else {
    alpha <- par[1]
    beta  <- par[2:(1 + k)]
    phi   <- phi_fixed
    sigma <- exp(par[2 + k])
    s     <- exp(par[3 + k])
  }
  
  # FTB contribution
  ll_ftb <- 0
  if (length(base_ctx$ftb_idx) > 0L) {
    B_ftb <- base_ctx$B0_full[base_ctx$ftb_idx, , drop = FALSE]
    mu_ftb <- alpha + drop(B_ftb %*% beta)
    ll_ftb <- sum(ald_logpdf(y[base_ctx$ftb_idx], mu_ftb, tau, s))
  }
  
  # sPTB contribution by deterministic quadrature
  ll_ptb <- 0
  if (length(base_ctx$ptb_idx) > 0L) {
    t_ptb <- base_ctx$t[base_ctx$ptb_idx]
    y_ptb <- y[base_ctx$ptb_idx]
    
    subj_by_node <- matrix(NA_real_,
                           nrow = length(base_ctx$ptb_ids),
                           ncol = quad_rule$nq)
    
    for (m in seq_len(quad_rule$nq)) {
      shift_m <- phi + sigma * quad_rule$nodes[m]
      B_ptb_m <- make_B(t_ptb + shift_m, base_ctx$spline_spec)
      mu_ptb_m <- alpha + drop(B_ptb_m %*% beta)
      
      ll_row_m <- ald_logpdf(y_ptb, mu_ptb_m, tau, s)
      subj_sum_m <- rowsum(ll_row_m, base_ctx$ptb_subject_factor, reorder = FALSE)
      
      subj_by_node[, m] <- as.numeric(subj_sum_m) + log(quad_rule$weights[m])
    }
    
    ll_ptb <- sum(apply(subj_by_node, 1, log_sum_exp))
  }
  
  nll <- -(ll_ftb + ll_ptb)
  if (!is.finite(nll)) nll <- 1e100
  nll
}

fit_hetero_mle <- function(base_ctx, y, tau, quad_rule,
                           start_common = NULL,
                           start_sigma_grid = SIGMA_INIT_GRID,
                           phi_fixed = NULL,
                           maxit = OPT_MAXIT_HET) {
  k <- base_ctx$n_basis
  
  if (is.null(start_common)) {
    start_common <- fit_common_mle(base_ctx, y, tau, phi_fixed = phi_fixed)
    if (is.null(start_common)) return(NULL)
  }
  
  cand <- lapply(start_sigma_grid, function(sig0) {
    if (is.null(phi_fixed)) {
      par0 <- c(start_common$alpha, start_common$beta,
                start_common$phi, log(sig0), log(start_common$s))
    } else {
      par0 <- c(start_common$alpha, start_common$beta,
                log(sig0), log(start_common$s))
    }
    val0 <- hetero_negloglik(par0, base_ctx, y, tau, quad_rule, phi_fixed = phi_fixed)
    list(par0 = par0, val0 = val0)
  })
  vals <- sapply(cand, `[[`, "val0")
  if (!any(is.finite(vals))) return(NULL)
  par0 <- cand[[which.min(vals)]]$par0
  
  if (is.null(phi_fixed)) {
    lower <- c(rep(-Inf, 1 + k), PHI_LOWER, log(SIGMA_LOWER), log(SCALE_LOWER))
    upper <- c(rep( Inf, 1 + k), PHI_UPPER, log(SIGMA_UPPER), log(SCALE_UPPER))
  } else {
    lower <- c(rep(-Inf, 1 + k), log(SIGMA_LOWER), log(SCALE_LOWER))
    upper <- c(rep( Inf, 1 + k), log(SIGMA_UPPER), log(SCALE_UPPER))
  }
  
  opt <- tryCatch(
    optim(
      par = par0,
      fn = hetero_negloglik,
      base_ctx = base_ctx,
      y = y,
      tau = tau,
      quad_rule = quad_rule,
      phi_fixed = phi_fixed,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = maxit)
    ),
    error = function(e) NULL
  )
  if (is.null(opt) || !is.finite(opt$value)) return(NULL)
  
  if (is.null(phi_fixed)) {
    alpha <- opt$par[1]
    beta  <- opt$par[2:(1 + base_ctx$n_basis)]
    phi   <- opt$par[2 + base_ctx$n_basis]
    sigma <- exp(opt$par[3 + base_ctx$n_basis])
    s     <- exp(opt$par[4 + base_ctx$n_basis])
  } else {
    alpha <- opt$par[1]
    beta  <- opt$par[2:(1 + base_ctx$n_basis)]
    phi   <- phi_fixed
    sigma <- exp(opt$par[2 + base_ctx$n_basis])
    s     <- exp(opt$par[3 + base_ctx$n_basis])
  }
  
  loglik <- -hetero_negloglik(opt$par, base_ctx, y, tau, quad_rule, phi_fixed = phi_fixed)
  
  list(
    alpha = alpha,
    beta = beta,
    phi = phi,
    sigma = sigma,
    sigma2 = sigma^2,
    s = s,
    loglik = loglik,
    conv = opt$convergence
  )
}

posterior_phi_i_eb <- function(fit_het, base_ctx, y, tau, quad_rule) {
  if (length(base_ctx$ptb_idx) == 0L) return(setNames(numeric(0), character(0)))
  
  alpha <- fit_het$alpha
  beta  <- fit_het$beta
  phi   <- fit_het$phi
  sigma <- fit_het$sigma
  s     <- fit_het$s
  
  t_ptb <- base_ctx$t[base_ctx$ptb_idx]
  y_ptb <- y[base_ctx$ptb_idx]
  
  subj_by_node <- matrix(NA_real_,
                         nrow = length(base_ctx$ptb_ids),
                         ncol = quad_rule$nq)
  
  for (m in seq_len(quad_rule$nq)) {
    shift_m <- phi + sigma * quad_rule$nodes[m]
    B_ptb_m <- make_B(t_ptb + shift_m, base_ctx$spline_spec)
    mu_ptb_m <- alpha + drop(B_ptb_m %*% beta)
    
    ll_row_m <- ald_logpdf(y_ptb, mu_ptb_m, tau, s)
    subj_sum_m <- rowsum(ll_row_m, base_ctx$ptb_subject_factor, reorder = FALSE)
    
    subj_by_node[, m] <- as.numeric(subj_sum_m) + log(quad_rule$weights[m])
  }
  
  eb <- numeric(length(base_ctx$ptb_ids))
  for (i in seq_len(length(base_ctx$ptb_ids))) {
    z <- subj_by_node[i, ]
    z <- z - log_sum_exp(z)
    w <- exp(z)
    eb[i] <- phi + sigma * sum(w * quad_rule$nodes)
  }
  setNames(eb, base_ctx$ptb_ids)
}

# -----------------------------
# Simulation from fitted models
# -----------------------------
simulate_from_common_fit <- function(base_ctx, fit_common, tau, phi_override = NULL) {
  phi_use <- if (is.null(phi_override)) fit_common$phi else phi_override
  B <- make_B(base_ctx$t + phi_use * base_ctx$g, base_ctx$spline_spec)
  mu <- fit_common$alpha + drop(B %*% fit_common$beta)
  mu + rALD(length(mu), tau = tau, s = fit_common$s)
}

simulate_from_hetero_fit <- function(base_ctx, fit_het, tau, dist_spec) {
  y_star <- numeric(length(base_ctx$t))
  
  # FTB
  if (length(base_ctx$ftb_idx) > 0L) {
    B_ftb <- base_ctx$B0_full[base_ctx$ftb_idx, , drop = FALSE]
    mu_ftb <- fit_het$alpha + drop(B_ftb %*% fit_het$beta)
    y_star[base_ctx$ftb_idx] <- mu_ftb + rALD(length(mu_ftb), tau = tau, s = fit_het$s)
  }
  
  # sPTB
  if (length(base_ctx$ptb_ids) > 0L) {
    u_i <- r_random_effect_std(length(base_ctx$ptb_ids), dist_spec)
    phi_i <- fit_het$phi + fit_het$sigma * u_i
    names(phi_i) <- base_ctx$ptb_ids
    
    for (sid in base_ctx$ptb_ids) {
      idx <- base_ctx$subj_index[[sid]]
      B_i <- make_B(base_ctx$t[idx] + phi_i[sid], base_ctx$spline_spec)
      mu_i <- fit_het$alpha + drop(B_i %*% fit_het$beta)
      y_star[idx] <- mu_i + rALD(length(idx), tau = tau, s = fit_het$s)
    }
  }
  
  y_star
}

# -----------------------------
# Sensitivity analysis for one biomarker/tau
# -----------------------------
analyze_one_task <- function(dat, yvar, tau,
                             dist_specs = DIST_SPECS,
                             id_col = "Participant_ID",
                             t_col = "GA",
                             g_col = "G",
                             df_spline = DF_SPLINE,
                             phi_grid = PHI_GRID,
                             nq_dist = NQ_DIST,
                             B_shift_common_test = B_SHIFT_COMMON_TEST,
                             B_shift_common_ci = B_SHIFT_COMMON_CI,
                             B_shift_het_test = B_SHIFT_HET_TEST,
                             B_shift_het_ci = B_SHIFT_HET_CI,
                             B_het_test = B_HET_TEST,
                             B_sig2_ci = B_SIG2_CI,
                             ncores_task = NCORES_TASK,
                             seed = BASE_SEED) {
  y <- dat[[yvar]]
  
  # Fixed spline specification reused throughout this task
  lower_extra <- min(phi_grid, na.rm = TRUE)
  upper_extra <- max(phi_grid, na.rm = TRUE)
  spline_spec <- build_spline_spec(
    x_base = dat[[t_col]],
    df_spline = df_spline,
    lower_extra = lower_extra,
    upper_extra = upper_extra
  )
  
  base_ctx <- build_base_context(
    dat = dat,
    yvar = yvar,
    spline_spec = spline_spec,
    phi_grid = phi_grid,
    id_col = id_col,
    t_col = t_col,
    g_col = g_col
  )
  
  n_obs <- nrow(dat)
  n_id  <- length(unique(dat[[id_col]]))
  n_ptb <- length(unique(dat[[id_col]][dat[[g_col]] == 1L]))
  n_ftb <- length(unique(dat[[id_col]][dat[[g_col]] == 0L]))
  
  # ============================================================
  # Reference common-shift model (distribution-free across specs)
  # ============================================================
  fit_common_alt <- fit_common_mle(base_ctx, y, tau, phi_fixed = NULL)
  if (is.null(fit_common_alt)) stop("Common-shift alternative fit failed.")
  
  fit_common_null_shift <- fit_common_mle(
    base_ctx, y, tau, phi_fixed = 0, start_fit = fit_common_alt
  )
  if (is.null(fit_common_null_shift)) stop("Common-shift null (phi=0) fit failed.")
  
  ll_common_alt <- fit_common_alt$loglik
  ll_common_null_shift <- fit_common_null_shift$loglik
  LR_shift_common_obs <- max(0, 2 * (ll_common_alt - ll_common_null_shift))
  
  # Parametric bootstrap p-value for common-model H0_shift
  common_shift_test_seeds <- seed + 100000L + seq_len(B_shift_common_test)
  common_shift_test_fun <- function(b) {
    set.seed(common_shift_test_seeds[b])
    y_star <- simulate_from_common_fit(base_ctx, fit_common_null_shift, tau, phi_override = 0)
    
    fit0_b <- fit_common_mle(base_ctx, y_star, tau, phi_fixed = 0, start_fit = fit_common_null_shift)
    if (is.null(fit0_b)) return(NA_real_)
    
    fit1_b <- fit_common_mle(base_ctx, y_star, tau, phi_fixed = NULL, start_fit = fit_common_alt)
    if (is.null(fit1_b)) return(NA_real_)
    
    max(0, 2 * (fit1_b$loglik - fit0_b$loglik))
  }
  
  LR_shift_common_boot <- as.numeric(unlist(parallel_lapply_safe(
    X = seq_len(B_shift_common_test),
    FUN = common_shift_test_fun,
    mc.cores = ncores_task
  )))
  p_shift_common <- mean(LR_shift_common_boot >= LR_shift_common_obs, na.rm = TRUE)
  
  # Parametric bootstrap CI for common-model phi
  common_shift_ci_seeds <- seed + 200000L + seq_len(B_shift_common_ci)
  common_shift_ci_fun <- function(b) {
    set.seed(common_shift_ci_seeds[b])
    y_star <- simulate_from_common_fit(base_ctx, fit_common_alt, tau)
    fit_b <- fit_common_mle(base_ctx, y_star, tau, phi_fixed = NULL, start_fit = fit_common_alt)
    if (is.null(fit_b)) return(NA_real_)
    fit_b$phi
  }
  
  phi_common_boot <- as.numeric(unlist(parallel_lapply_safe(
    X = seq_len(B_shift_common_ci),
    FUN = common_shift_ci_fun,
    mc.cores = ncores_task
  )))
  phi_common_ci <- safe_quantile(phi_common_boot, c(0.025, 0.975))
  
  # ============================================================
  # Sensitivity loop over random-effect distributions
  # ============================================================
  out_rows <- vector("list", length(dist_specs))
  
  for (d in seq_along(dist_specs)) {
    dist_spec <- dist_specs[[d]]
    quad_rule <- build_quad_rule(dist_spec, nq = nq_dist)
    
    message("  Fitting random-effect distribution: ", dist_spec$name)
    
    row_out <- list(
      biomarker = yvar,
      tau = tau,
      dist_name = dist_spec$name,
      dist_family = dist_spec$family,
      dist_param = dist_spec$param,
      n_obs = n_obs,
      n_id = n_id,
      n_ptb = n_ptb,
      n_ftb = n_ftb,
      
      # common-model reference
      phi_common_hat = fit_common_alt$phi,
      phi_common_ci_low = phi_common_ci[1],
      phi_common_ci_high = phi_common_ci[2],
      ll_common_alt = ll_common_alt,
      ll_common_null_shift = ll_common_null_shift,
      LR_shift_common_obs = LR_shift_common_obs,
      p_shift_common = p_shift_common,
      s_common_hat = fit_common_alt$s,
      
      # spec-specific hetero outputs initialized
      phi_het_hat = NA_real_,
      phi_het_ci_low = NA_real_,
      phi_het_ci_high = NA_real_,
      sigma_hat = NA_real_,
      sig2_hat = NA_real_,
      sig2_ci_low = NA_real_,
      sig2_ci_high = NA_real_,
      ll_het_alt = NA_real_,
      ll_het_null_shift = NA_real_,
      LR_shift_het_obs = NA_real_,
      p_shift_het = NA_real_,
      LR_het_obs = NA_real_,
      p_het = NA_real_,
      s_het_hat = NA_real_,
      eb_phi_var_hat = NA_real_,
      n_nonzero_dev_phi_eb = NA_integer_,
      max_abs_dev_phi_eb = NA_real_,
      ncores_task = ncores_task,
      error = NA_character_
    )
    
    if (n_ptb <= 1L) {
      row_out$error <- "Too few sPTB subjects for heterogeneous model."
      out_rows[[d]] <- row_out
      next
    }
    
    fit_het_alt <- tryCatch(
      fit_hetero_mle(
        base_ctx = base_ctx,
        y = y,
        tau = tau,
        quad_rule = quad_rule,
        start_common = fit_common_alt,
        phi_fixed = NULL
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit_het_alt)) {
      row_out$error <- "Heterogeneous alternative fit failed."
      out_rows[[d]] <- row_out
      next
    }
    
    fit_het_null_shift <- tryCatch(
      fit_hetero_mle(
        base_ctx = base_ctx,
        y = y,
        tau = tau,
        quad_rule = quad_rule,
        start_common = fit_common_null_shift,
        phi_fixed = 0
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit_het_null_shift)) {
      row_out$error <- "Heterogeneous null-shift fit failed."
      out_rows[[d]] <- row_out
      next
    }
    
    # Point estimates
    row_out$phi_het_hat <- fit_het_alt$phi
    row_out$sigma_hat <- fit_het_alt$sigma
    row_out$sig2_hat <- fit_het_alt$sigma2
    row_out$ll_het_alt <- fit_het_alt$loglik
    row_out$ll_het_null_shift <- fit_het_null_shift$loglik
    row_out$s_het_hat <- fit_het_alt$s
    row_out$LR_shift_het_obs <- max(0, 2 * (fit_het_alt$loglik - fit_het_null_shift$loglik))
    row_out$LR_het_obs <- max(0, 2 * (fit_het_alt$loglik - ll_common_alt))
    
    # EB diagnostics
    eb_phi <- posterior_phi_i_eb(fit_het_alt, base_ctx, y, tau, quad_rule)
    dev_eb <- eb_phi - fit_het_alt$phi
    row_out$eb_phi_var_hat <- if (length(eb_phi) > 1L) var(eb_phi) else 0
    row_out$n_nonzero_dev_phi_eb <- sum(abs(dev_eb) > 1e-8)
    row_out$max_abs_dev_phi_eb <- if (length(dev_eb) > 0L) max(abs(dev_eb)) else 0
    
    # ----------------------------------------------------------
    # Spec-specific bootstrap test for H0_shift: phi = 0
    # ----------------------------------------------------------
    het_shift_test_seeds <- seed + 300000L + d * 10000L + seq_len(B_shift_het_test)
    het_shift_test_fun <- function(b) {
      set.seed(het_shift_test_seeds[b])
      
      y_star <- simulate_from_hetero_fit(
        base_ctx = base_ctx,
        fit_het = fit_het_null_shift,
        tau = tau,
        dist_spec = dist_spec
      )
      
      fit0_b <- fit_hetero_mle(
        base_ctx = base_ctx,
        y = y_star,
        tau = tau,
        quad_rule = quad_rule,
        start_common = fit_common_null_shift,
        phi_fixed = 0
      )
      if (is.null(fit0_b)) return(NA_real_)
      
      fit1_b <- fit_hetero_mle(
        base_ctx = base_ctx,
        y = y_star,
        tau = tau,
        quad_rule = quad_rule,
        start_common = fit_common_alt,
        phi_fixed = NULL
      )
      if (is.null(fit1_b)) return(NA_real_)
      
      max(0, 2 * (fit1_b$loglik - fit0_b$loglik))
    }
    
    LR_shift_het_boot <- as.numeric(unlist(parallel_lapply_safe(
      X = seq_len(B_shift_het_test),
      FUN = het_shift_test_fun,
      mc.cores = ncores_task
    )))
    row_out$p_shift_het <- mean(LR_shift_het_boot >= row_out$LR_shift_het_obs, na.rm = TRUE)
    
    # ----------------------------------------------------------
    # Spec-specific bootstrap CI for phi under hetero model
    # ----------------------------------------------------------
    het_phi_ci_seeds <- seed + 400000L + d * 10000L + seq_len(B_shift_het_ci)
    het_phi_ci_fun <- function(b) {
      set.seed(het_phi_ci_seeds[b])
      
      y_star <- simulate_from_hetero_fit(
        base_ctx = base_ctx,
        fit_het = fit_het_alt,
        tau = tau,
        dist_spec = dist_spec
      )
      
      fit_b <- fit_hetero_mle(
        base_ctx = base_ctx,
        y = y_star,
        tau = tau,
        quad_rule = quad_rule,
        start_common = fit_common_alt,
        phi_fixed = NULL
      )
      if (is.null(fit_b)) return(NA_real_)
      fit_b$phi
    }
    
    phi_het_boot <- as.numeric(unlist(parallel_lapply_safe(
      X = seq_len(B_shift_het_ci),
      FUN = het_phi_ci_fun,
      mc.cores = ncores_task
    )))
    phi_het_ci <- safe_quantile(phi_het_boot, c(0.025, 0.975))
    row_out$phi_het_ci_low <- phi_het_ci[1]
    row_out$phi_het_ci_high <- phi_het_ci[2]
    
    # ----------------------------------------------------------
    # Spec-specific bootstrap test for H0_het: sigma^2 = 0
    # ----------------------------------------------------------
    het_test_seeds <- seed + 500000L + d * 10000L + seq_len(B_het_test)
    het_test_fun <- function(b) {
      set.seed(het_test_seeds[b])
      
      y_star <- simulate_from_common_fit(
        base_ctx = base_ctx,
        fit_common = fit_common_alt,
        tau = tau
      )
      
      fit0_b <- fit_common_mle(
        base_ctx = base_ctx,
        y = y_star,
        tau = tau,
        phi_fixed = NULL,
        start_fit = fit_common_alt
      )
      if (is.null(fit0_b)) return(NA_real_)
      
      fit1_b <- fit_hetero_mle(
        base_ctx = base_ctx,
        y = y_star,
        tau = tau,
        quad_rule = quad_rule,
        start_common = fit0_b,
        phi_fixed = NULL
      )
      if (is.null(fit1_b)) return(NA_real_)
      
      max(0, 2 * (fit1_b$loglik - fit0_b$loglik))
    }
    
    LR_het_boot <- as.numeric(unlist(parallel_lapply_safe(
      X = seq_len(B_het_test),
      FUN = het_test_fun,
      mc.cores = ncores_task
    )))
    row_out$p_het <- mean(LR_het_boot >= row_out$LR_het_obs, na.rm = TRUE)
    
    # ----------------------------------------------------------
    # Spec-specific bootstrap CI for sigma^2
    # ----------------------------------------------------------
    sig2_ci_seeds <- seed + 600000L + d * 10000L + seq_len(B_SIG2_CI)
    sig2_ci_fun <- function(b) {
      set.seed(sig2_ci_seeds[b])
      
      y_star <- simulate_from_hetero_fit(
        base_ctx = base_ctx,
        fit_het = fit_het_alt,
        tau = tau,
        dist_spec = dist_spec
      )
      
      fit0_b <- fit_common_mle(
        base_ctx = base_ctx,
        y = y_star,
        tau = tau,
        phi_fixed = NULL,
        start_fit = fit_common_alt
      )
      if (is.null(fit0_b)) return(NA_real_)
      
      fit1_b <- fit_hetero_mle(
        base_ctx = base_ctx,
        y = y_star,
        tau = tau,
        quad_rule = quad_rule,
        start_common = fit0_b,
        phi_fixed = NULL
      )
      if (is.null(fit1_b)) return(NA_real_)
      
      fit1_b$sigma2
    }
    
    sig2_boot <- as.numeric(unlist(parallel_lapply_safe(
      X = seq_len(B_SIG2_CI),
      FUN = sig2_ci_fun,
      mc.cores = ncores_task
    )))
    sig2_ci <- safe_quantile(sig2_boot, c(0.025, 0.975))
    row_out$sig2_ci_low <- sig2_ci[1]
    row_out$sig2_ci_high <- sig2_ci[2]
    
    out_rows[[d]] <- row_out
  }
  
  bind_rows(out_rows)
}

# -----------------------------
# Run one task
# -----------------------------
res <- tryCatch(
  analyze_one_task(
    dat = df,
    yvar = yvar,
    tau = tau,
    dist_specs = DIST_SPECS,
    df_spline = DF_SPLINE,
    phi_grid = PHI_GRID,
    nq_dist = NQ_DIST,
    B_shift_common_test = B_SHIFT_COMMON_TEST,
    B_shift_common_ci = B_SHIFT_COMMON_CI,
    B_shift_het_test = B_SHIFT_HET_TEST,
    B_shift_het_ci = B_SHIFT_HET_CI,
    B_het_test = B_HET_TEST,
    B_sig2_ci = B_SIG2_CI,
    ncores_task = NCORES_TASK,
    seed = BASE_SEED + task_id
  ),
  error = function(e) {
    data.frame(
      biomarker = yvar,
      tau = tau,
      dist_name = NA_character_,
      error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  }
)

write.csv(res, part_file, row.names = FALSE)
message("Wrote: ", part_file)