#!/usr/bin/env Rscript

# Multi-response analyses for tex/multi_response_three_analyses.tex.
#
# The script performs exactly the three analyses in that note:
#   1) multi-response gestational-age quantile spline tests;
#   2) multi-response forward-time common phase-shift analysis;
#   3) multi-response FTB-only phase-variance analysis.
#
# Outputs:
#   Test_results_multi_response_three_analyses/*.csv
#   Test_results_multi_response_three_analyses/*.rds
#   tex/multi_response_three_analyses_results.tex

rm(list = ls())

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(quantreg)
  library(splines)
})

# ---------------------------------------------------------------------------
# Paths and settings
# ---------------------------------------------------------------------------

script_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
  error = function(e) NA_character_
)

if (!is.na(script_path)) {
  project_root <- normalizePath(file.path(dirname(script_path), ".."),
                                winslash = "/", mustWork = TRUE)
} else {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  project_root <- if (basename(wd) == "codes") dirname(wd) else wd
}

DATA_PATH <- Sys.getenv(
  "MR_DATA_PATH",
  unset = file.path(project_root, "data", "V1_V7_longformat_k1.csv")
)
OUT_DIR <- Sys.getenv(
  "MR_OUT_DIR",
  unset = file.path(project_root, "Test_results_multi_response_three_analyses")
)
TEX_OUT <- Sys.getenv(
  "MR_TEX_OUT",
  unset = file.path(project_root, "tex", "multi_response_three_analyses_results.tex")
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(TEX_OUT), showWarnings = FALSE, recursive = TRUE)

BIOMARKERS <- c("SWS", "Midband", "Intercept", "AC", "Slope")
TAUS <- as.numeric(strsplit(Sys.getenv("MR_TAUS", "0.20,0.50,0.80"), ",")[[1]])
DF_SPLINE <- as.integer(Sys.getenv("MR_DF_SPLINE", "3"))
B_BOOT <- as.integer(Sys.getenv("MR_B_BOOT", "300"))
SEED0 <- as.integer(Sys.getenv("MR_SEED", "20260515"))

PHI_GRID <- seq(
  as.numeric(Sys.getenv("MR_PHI_LOWER", "-6")),
  as.numeric(Sys.getenv("MR_PHI_UPPER", "6")),
  by = as.numeric(Sys.getenv("MR_PHI_STEP", "0.25"))
)
PHI_REFINE_STEP <- as.numeric(Sys.getenv("MR_PHI_REFINE_STEP", "0.05"))
PHI_REFINE_RADIUS <- as.numeric(Sys.getenv("MR_PHI_REFINE_RADIUS", "0.50"))
H_SPARSE <- as.numeric(Sys.getenv("MR_H_SPARSE", "0.05"))
DERIV_DELTA <- as.numeric(Sys.getenv("MR_DERIV_DELTA", "0.05"))
MIN_SPARSITY <- as.numeric(Sys.getenv("MR_MIN_SPARSITY", "1e-4"))

set.seed(SEED0)

# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

fmt_num <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(
    is.na(x),
    "--",
    ifelse(abs(x) < 0.5 * 10^(-digits),
           paste0("0.", paste(rep("0", digits), collapse = "")),
           formatC(x, format = "f", digits = digits))
  )
}

fmt_p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(
    is.na(x),
    "--",
    ifelse(x < 0.001, "<0.001", formatC(x, format = "f", digits = 3))
  )
}

tex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("&", "\\\\&", x)
  x
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

rho_tau <- function(u, tau) u * (tau - (u < 0))

safe_mad <- function(x) {
  x <- as_num(x)
  med <- median(x, na.rm = TRUE)
  s <- mad(x, center = med, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(s) || s <= 1e-10) {
    qs <- quantile(x, probs = c(0.25, 0.75), names = FALSE, na.rm = TRUE)
    s <- diff(qs) / 1.349
  }
  if (!is.finite(s) || s <= 1e-10) s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 1e-10) s <- 1
  s
}

ga_multiplier <- function(x) {
  mx <- suppressWarnings(max(as_num(x), na.rm = TRUE))
  if (is.finite(mx) && mx <= 2) 39 else 1
}

outcome_to_g <- function(dat) {
  raw <- as.character(dat$Outcome)
  g <- ifelse(raw %in% c("1", "PTB", "sPTB", "ptb", "sptb"), 1L,
              ifelse(raw %in% c("0", "FTB", "ftb"), 0L, NA_integer_))
  if (all(!is.na(g)) && length(unique(g)) == 2) return(g)

  if ("mPTB" %in% names(dat)) {
    raw2 <- as.character(dat$mPTB)
    g2 <- ifelse(raw2 %in% c("1", "PTB", "sPTB", "ptb", "sptb"), 1L,
                 ifelse(raw2 %in% c("0", "FTB", "ftb"), 0L, NA_integer_))
    if (sum(!is.na(g2)) > 0 && length(unique(g2[!is.na(g2)])) == 2) return(g2)
  }
  stop("Could not infer binary PTB/FTB group from Outcome or mPTB.")
}

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

if (!file.exists(DATA_PATH)) stop("Missing data file: ", DATA_PATH)

raw <- read_csv(DATA_PATH, show_col_types = FALSE)
ga_mult <- ga_multiplier(raw$GA)

dat <- raw %>%
  mutate(
    Participant_ID = as.character(Participant_ID),
    Visit_ID = if ("Visit_ID" %in% names(.)) as.character(Visit_ID) else as.character(row_number()),
    G = outcome_to_g(.),
    GA = as_num(GA) * ga_mult,
    GA_delivery = as_num(GA_delivery) * ga_mult
  ) %>%
  filter(!is.na(Participant_ID), !is.na(G), !is.na(GA))

for (v in BIOMARKERS) dat[[v]] <- as_num(dat[[v]])

if (any(tapply(dat$G, dat$Participant_ID, function(z) length(unique(z)) > 1))) {
  stop("G is not constant within Participant_ID.")
}

dat <- dat %>% arrange(Participant_ID, GA)
dat_ftb <- dat %>% filter(G == 0)
dat_ptb <- dat %>% filter(G == 1)

data_summary <- tibble(
  quantity = c("Subjects", "FTB subjects", "PTB subjects",
               "Visits", "FTB visits", "PTB visits"),
  value = c(
    n_distinct(dat$Participant_ID),
    n_distinct(dat_ftb$Participant_ID),
    n_distinct(dat_ptb$Participant_ID),
    nrow(dat),
    nrow(dat_ftb),
    nrow(dat_ptb)
  )
)
write_csv(data_summary, file.path(OUT_DIR, "data_summary.csv"))

message("Project root: ", project_root)
message("Rows: ", nrow(dat), "; subjects: ", n_distinct(dat$Participant_ID))
message("Bootstrap replicates: ", B_BOOT)

# ---------------------------------------------------------------------------
# Analysis 1. Multi-response GA spline QR
# ---------------------------------------------------------------------------

select_group_terms <- function(coef_names) {
  coef_names[
    coef_names == "G" |
      grepl("^G:ns\\(GA", coef_names) |
      (grepl(":G$", coef_names) & grepl("ns\\(GA", coef_names))
  ]
}

make_param_name <- function(response, tau, term) {
  paste(response, sprintf("tau%.2f", tau), term, sep = "__")
}

fit_ga_group_vector <- function(d, taus, biomarkers, param_specs = NULL) {
  beta <- numeric(0)
  spec_list <- list()

  for (tau in taus) {
    for (v in biomarkers) {
      needed <- c("GA", "G", "Participant_ID", v)
      ds <- d[complete.cases(d[, needed, drop = FALSE]), , drop = FALSE]

      if (!is.null(param_specs)) {
        target_terms <- param_specs$term[param_specs$response == v &
                                           abs(param_specs$tau - tau) < 1e-12]
      } else {
        target_terms <- NULL
      }

      if (nrow(ds) < 30 || length(unique(ds$G)) < 2) {
        if (length(target_terms) > 0) {
          vals <- rep(NA_real_, length(target_terms))
          names(vals) <- make_param_name(v, tau, target_terms)
          beta <- c(beta, vals)
        }
        next
      }

      form <- as.formula(sprintf("%s ~ G * ns(GA, df = %d)", v, DF_SPLINE))
      fit <- tryCatch(
        suppressWarnings(rq(form, tau = tau, data = ds)),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        if (length(target_terms) > 0) {
          vals <- rep(NA_real_, length(target_terms))
          names(vals) <- make_param_name(v, tau, target_terms)
          beta <- c(beta, vals)
        }
        next
      }

      cf <- coef(fit)
      if (is.null(target_terms)) target_terms <- select_group_terms(names(cf))
      vals <- unname(cf[target_terms])
      vals[is.na(vals)] <- NA_real_
      names(vals) <- make_param_name(v, tau, target_terms)
      beta <- c(beta, vals)

      if (is.null(param_specs)) {
        spec_list[[length(spec_list) + 1]] <- tibble(
          tau = tau,
          response = v,
          term = target_terms,
          param_name = make_param_name(v, tau, target_terms)
        )
      }
    }
  }

  specs <- if (is.null(param_specs)) bind_rows(spec_list) else param_specs
  if (!is.null(param_specs)) beta <- beta[specs$param_name]
  list(beta = beta, specs = specs)
}

regularized_wald <- function(beta_hat, boot_mat, idx, tau_label) {
  idx <- as.integer(idx)
  b0 <- beta_hat[idx]
  M <- boot_mat[, idx, drop = FALSE]
  ok <- complete.cases(M)
  M <- M[ok, , drop = FALSE]
  p <- length(idx)

  if (nrow(M) < max(25, p + 5)) {
    return(tibble(
      tau = tau_label,
      n_biomarkers = length(BIOMARKERS),
      group_df = p,
      W = NA_real_,
      p_boot = NA_real_,
      p_chisq = NA_real_,
      eig_rank = NA_integer_,
      B_ok = nrow(M),
      note = "Too few complete bootstrap replicates"
    ))
  }

  centered <- sweep(M, 2, b0, "-")
  S <- cov(centered)
  S <- (S + t(S)) / 2
  d <- diag(S)
  ridge <- median(d[is.finite(d) & d > 0], na.rm = TRUE)
  if (!is.finite(ridge) || ridge <= 0) ridge <- 1
  S <- S + diag(ridge * 1e-5, p)

  eg <- eigen(S, symmetric = TRUE)
  max_ev <- max(eg$values, na.rm = TRUE)
  keep <- eg$values > max(max_ev * 1e-8, 1e-12)
  rank <- sum(keep)
  Vinv <- eg$vectors[, keep, drop = FALSE] %*%
    diag(1 / eg$values[keep], nrow = rank) %*%
    t(eg$vectors[, keep, drop = FALSE])

  W <- as.numeric(t(b0) %*% Vinv %*% b0)
  W_star <- apply(centered, 1, function(z) as.numeric(t(z) %*% Vinv %*% z))
  W_star <- W_star[is.finite(W_star)]

  tibble(
    tau = tau_label,
    n_biomarkers = length(BIOMARKERS),
    group_df = p,
    W = W,
    p_boot = (1 + sum(W_star >= W)) / (1 + length(W_star)),
    p_chisq = pchisq(W, df = rank, lower.tail = FALSE),
    eig_rank = rank,
    B_ok = length(W_star),
    note = ""
  )
}

run_mr_ga <- function(d) {
  message("Analysis 1: observed multi-response GA spline fits")
  obs <- fit_ga_group_vector(d, TAUS, BIOMARKERS)
  beta_hat <- obs$beta
  specs <- obs$specs

  if (length(beta_hat) == 0 || anyNA(beta_hat)) {
    stop("Observed GA stacked fit failed.")
  }

  ids <- unique(d$Participant_ID)
  rows_by_id <- split(seq_len(nrow(d)), d$Participant_ID)
  boot_mat <- matrix(
    NA_real_,
    nrow = B_BOOT,
    ncol = length(beta_hat),
    dimnames = list(NULL, names(beta_hat))
  )

  set.seed(SEED0 + 101)
  message("Analysis 1: subject-level bootstrap")
  for (b in seq_len(B_BOOT)) {
    sid <- sample(ids, size = length(ids), replace = TRUE)
    db <- bind_rows(lapply(seq_along(sid), function(k) {
      z <- d[rows_by_id[[sid[k]]], , drop = FALSE]
      z$Participant_ID <- paste0(sid[k], "__boot__", k)
      z
    }))

    fit_b <- tryCatch(
      fit_ga_group_vector(db, TAUS, BIOMARKERS, param_specs = specs),
      error = function(e) NULL
    )
    if (!is.null(fit_b)) boot_mat[b, ] <- fit_b$beta[names(beta_hat)]
    if (b %% 50 == 0 || b == B_BOOT) {
      message("  GA bootstrap ", b, "/", B_BOOT)
    }
  }

  tests <- bind_rows(lapply(TAUS, function(tt) {
    idx <- which(abs(specs$tau - tt) < 1e-12)
    regularized_wald(beta_hat, boot_mat, idx, sprintf("%.2f", tt))
  }))
  tests <- bind_rows(
    tests,
    regularized_wald(beta_hat, boot_mat, seq_along(beta_hat), "Joint all quantiles")
  )

  term_se <- apply(boot_mat, 2, function(z) sd(z - beta_hat[names(z)], na.rm = TRUE))
  terms <- specs %>%
    mutate(
      estimate = as.numeric(beta_hat[param_name]),
      boot_se = as.numeric(term_se[param_name]),
      z = estimate / boot_se,
      p_normal = 2 * pnorm(abs(z), lower.tail = FALSE)
    )

  saveRDS(
    list(beta_hat = beta_hat, boot_mat = boot_mat, specs = specs,
         tests = tests, terms = terms),
    file.path(OUT_DIR, "analysis1_mr_ga_spline_bootstrap.rds")
  )
  write_csv(tests, file.path(OUT_DIR, "analysis1_mr_ga_spline_tests.csv"))
  write_csv(terms, file.path(OUT_DIR, "analysis1_mr_ga_spline_terms.csv"))
  tests
}

# ---------------------------------------------------------------------------
# Shared helpers for phase-shift analyses
# ---------------------------------------------------------------------------

standardization_from_ftb <- function(d) {
  centers <- sapply(d[, BIOMARKERS, drop = FALSE], median, na.rm = TRUE)
  scales <- sapply(d[, BIOMARKERS, drop = FALSE], safe_mad)
  list(centers = centers, scales = scales)
}

fit_std_reference_curves <- function(d_ftb, tau, std) {
  curves <- list()
  for (v in BIOMARKERS) {
    ds <- d_ftb %>%
      transmute(
        GA = GA,
        y = (.data[[v]] - std$centers[[v]]) / std$scales[[v]]
      ) %>%
      filter(is.finite(GA), is.finite(y))
    curves[[v]] <- tryCatch(
      suppressWarnings(rq(y ~ ns(GA, df = DF_SPLINE), tau = tau, data = ds)),
      error = function(e) NULL
    )
    if (is.null(curves[[v]])) stop("Failed reference curve for ", v, ", tau=", tau)
  }
  curves
}

pred_curve <- function(fit, t) {
  as.numeric(predict(fit, newdata = data.frame(GA = t)))
}

phase_loss <- function(d_target, curves, tau, phi, std) {
  total <- 0
  for (v in BIOMARKERS) {
    y <- (d_target[[v]] - std$centers[[v]]) / std$scales[[v]]
    keep <- is.finite(y) & is.finite(d_target$GA)
    if (!any(keep)) next
    q <- pred_curve(curves[[v]], d_target$GA[keep] + phi)
    total <- total + sum(rho_tau(y[keep] - q, tau), na.rm = TRUE)
  }
  total
}

grid_minimize <- function(fn, grid = PHI_GRID) {
  vals <- sapply(grid, fn)
  phi0 <- grid[which.min(vals)]
  lo <- max(min(grid), phi0 - PHI_REFINE_RADIUS)
  hi <- min(max(grid), phi0 + PHI_REFINE_RADIUS)
  fine <- seq(lo, hi, by = PHI_REFINE_STEP)
  vals2 <- sapply(fine, fn)
  list(phi = fine[which.min(vals2)], loss = min(vals2, na.rm = TRUE))
}

estimate_common_phi <- function(d_target, curves, tau, std) {
  grid_minimize(function(ph) phase_loss(d_target, curves, tau, ph, std))
}

estimate_subject_phis <- function(d_target, curves, tau, std, min_visits = 1L) {
  d_target <- d_target %>%
    filter(!is.na(Participant_ID), is.finite(GA))
  if (nrow(d_target) == 0) return(tibble())

  ids <- unique(d_target$Participant_ID)
  visit_counts <- tapply(d_target$Visit_ID, d_target$Participant_ID, function(z) length(unique(z)))
  keep_ids <- names(visit_counts)[visit_counts >= min_visits]
  d_target <- d_target[d_target$Participant_ID %in% keep_ids, , drop = FALSE]
  if (nrow(d_target) == 0) return(tibble())

  grid <- PHI_GRID
  loss_by_row <- matrix(0, nrow = nrow(d_target), ncol = length(grid))

  for (v in BIOMARKERS) {
    y <- (d_target[[v]] - std$centers[[v]]) / std$scales[[v]]
    ok <- is.finite(y)
    if (!any(ok)) next
    for (k in seq_along(grid)) {
      q <- pred_curve(curves[[v]], d_target$GA[ok] + grid[k])
      loss_by_row[ok, k] <- loss_by_row[ok, k] + rho_tau(y[ok] - q, tau)
    }
  }

  loss_by_subject <- rowsum(loss_by_row, group = d_target$Participant_ID, reorder = FALSE)
  best <- max.col(-loss_by_subject, ties.method = "first")
  tibble(
    Participant_ID = rownames(loss_by_subject),
    tau = tau,
    phi_i = grid[best],
    loss = loss_by_subject[cbind(seq_len(nrow(loss_by_subject)), best)],
    n_visits = as.integer(visit_counts[rownames(loss_by_subject)])
  )
}

sample_by_subject <- function(d, ids = unique(d$Participant_ID)) {
  rows_by_id <- split(seq_len(nrow(d)), d$Participant_ID)
  sid <- sample(ids, size = length(ids), replace = TRUE)
  bind_rows(lapply(seq_along(sid), function(k) {
    z <- d[rows_by_id[[sid[k]]], , drop = FALSE]
    z$Participant_ID <- paste0(sid[k], "__boot__", k)
    z
  }))
}

variance_boot_from_subject_shifts <- function(phi_i, B, seed) {
  phi_i <- phi_i[is.finite(phi_i)]
  if (length(phi_i) < 3) return(numeric(0))
  set.seed(seed)
  replicate(B, {
    z <- sample(phi_i, size = length(phi_i), replace = TRUE)
    var(z)
  })
}

# ---------------------------------------------------------------------------
# Analysis 2. Forward-time common phase shift
# ---------------------------------------------------------------------------

run_forward_phase <- function() {
  message("Analysis 2: multi-response forward-time common phase shift")
  rows <- list()
  subj_rows <- list()

  ids_ftb <- unique(dat_ftb$Participant_ID)
  ids_ptb <- unique(dat_ptb$Participant_ID)

  for (tau in TAUS) {
    std <- standardization_from_ftb(dat_ftb)
    curves <- fit_std_reference_curves(dat_ftb, tau, std)
    common <- estimate_common_phi(dat_ptb, curves, tau, std)
    subj <- estimate_subject_phis(dat_ptb, curves, tau, std, min_visits = 1L)
    subj_rows[[length(subj_rows) + 1]] <- subj
    sig2 <- var(subj$phi_i, na.rm = TRUE)

    phi_boot <- rep(NA_real_, B_BOOT)
    set.seed(SEED0 + 2000 + round(100 * tau))
    for (b in seq_len(B_BOOT)) {
      db_ftb <- sample_by_subject(dat_ftb, ids_ftb)
      db_ptb <- sample_by_subject(dat_ptb, ids_ptb)
      out <- tryCatch({
        std_b <- standardization_from_ftb(db_ftb)
        curves_b <- fit_std_reference_curves(db_ftb, tau, std_b)
        estimate_common_phi(db_ptb, curves_b, tau, std_b)$phi
      }, error = function(e) NA_real_)
      phi_boot[b] <- out
      if (b %% 50 == 0 || b == B_BOOT) {
        message("  phase tau=", tau, " bootstrap ", b, "/", B_BOOT)
      }
    }
    phi_boot <- phi_boot[is.finite(phi_boot)]
    sig2_boot <- variance_boot_from_subject_shifts(
      subj$phi_i,
      B_BOOT,
      SEED0 + 2100 + round(100 * tau)
    )

    ci_phi <- if (length(phi_boot) > 10) {
      as.numeric(quantile(phi_boot, c(0.025, 0.975), na.rm = TRUE))
    } else c(NA_real_, NA_real_)
    ci_sig2 <- if (length(sig2_boot) > 10) {
      as.numeric(quantile(sig2_boot, c(0.025, 0.975), na.rm = TRUE))
    } else c(NA_real_, NA_real_)
    p_shift <- if (length(phi_boot) > 10) {
      min(1, 2 * min(
        (1 + sum(phi_boot <= 0)) / (1 + length(phi_boot)),
        (1 + sum(phi_boot >= 0)) / (1 + length(phi_boot))
      ))
    } else NA_real_
    p_sig2_le0 <- if (length(sig2_boot) > 10) {
      (1 + sum(sig2_boot <= 1e-10)) / (1 + length(sig2_boot))
    } else NA_real_

    rows[[length(rows) + 1]] <- tibble(
      tau = tau,
      phi_hat_weeks = common$phi,
      phi_ci_low = ci_phi[1],
      phi_ci_high = ci_phi[2],
      phi_hat_days = 7 * common$phi,
      p_shift_boot = p_shift,
      sigma2_phi_hat = sig2,
      sigma2_phi_ci_low = ci_sig2[1],
      sigma2_phi_ci_high = ci_sig2[2],
      p_sigma2_le0_boot = p_sig2_le0,
      B_phi_ok = length(phi_boot),
      n_ptb_subjects_for_phi_i = nrow(subj),
      note = "Common shift from FTB reference curves; variance is descriptive subject-specific profile-shift variance."
    )
  }

  res <- bind_rows(rows)
  subjects <- bind_rows(subj_rows)
  write_csv(res, file.path(OUT_DIR, "analysis2_forward_time_phase_shift.csv"))
  write_csv(subjects, file.path(OUT_DIR, "analysis2_forward_time_phase_subject_shifts.csv"))
  res
}

# ---------------------------------------------------------------------------
# Analysis 3. FTB-only phase variance
# ---------------------------------------------------------------------------

fit_raw_quantile_curves <- function(d_ftb, tau) {
  curves <- list()
  for (v in BIOMARKERS) {
    ds <- d_ftb %>%
      select(GA, all_of(v)) %>%
      filter(is.finite(GA), is.finite(.data[[v]]))
    curves[[v]] <- list(
      q0 = tryCatch(suppressWarnings(rq(as.formula(paste0(v, " ~ ns(GA, df = ", DF_SPLINE, ")")),
                                      tau = tau, data = ds)), error = function(e) NULL),
      qm = tryCatch(suppressWarnings(rq(as.formula(paste0(v, " ~ ns(GA, df = ", DF_SPLINE, ")")),
                                      tau = tau - H_SPARSE, data = ds)), error = function(e) NULL),
      qp = tryCatch(suppressWarnings(rq(as.formula(paste0(v, " ~ ns(GA, df = ", DF_SPLINE, ")")),
                                      tau = tau + H_SPARSE, data = ds)), error = function(e) NULL)
    )
    if (any(sapply(curves[[v]], is.null))) stop("Failed raw curves for ", v, ", tau=", tau)
  }
  curves
}

deriv_curve <- function(fit, t) {
  (pred_curve(fit, t + DERIV_DELTA) - pred_curve(fit, t - DERIV_DELTA)) /
    (2 * DERIV_DELTA)
}

make_hit_lambda_long <- function(d_ftb, tau, curves) {
  rows <- list()
  for (v in BIOMARKERS) {
    keep <- is.finite(d_ftb[[v]]) & is.finite(d_ftb$GA)
    ds <- d_ftb[keep, , drop = FALSE]
    q0 <- pred_curve(curves[[v]]$q0, ds$GA)
    qm <- pred_curve(curves[[v]]$qm, ds$GA)
    qp <- pred_curve(curves[[v]]$qp, ds$GA)
    sparsity <- pmax(abs((qp - qm) / (2 * H_SPARSE)), MIN_SPARSITY)
    ghat <- 1 / sparsity
    dq <- deriv_curve(curves[[v]]$q0, ds$GA)

    rows[[length(rows) + 1]] <- tibble(
      Participant_ID = ds$Participant_ID,
      Visit_ID = ds$Visit_ID,
      biomarker = v,
      GA = ds$GA,
      H = as.numeric(ds[[v]] <= q0) - tau,
      lambda = ghat * dq
    )
  }
  bind_rows(rows)
}

pair_contrib_mr <- function(d) {
  if (nrow(d) < 2 || n_distinct(d$Visit_ID) < 2) {
    return(tibble(N = 0, D = 0, npair = 0))
  }
  cmb <- combn(seq_len(nrow(d)), 2)
  cross_time <- d$Visit_ID[cmb[1, ]] != d$Visit_ID[cmb[2, ]]
  ll <- d$lambda[cmb[1, ]] * d$lambda[cmb[2, ]]
  hh <- d$H[cmb[1, ]] * d$H[cmb[2, ]]
  keep <- cross_time & is.finite(ll) & is.finite(hh) & abs(ll) > 1e-10
  tibble(
    N = sum(ll[keep] * hh[keep]),
    D = sum(ll[keep]^2),
    npair = sum(keep)
  )
}

estimate_mr_hit_moment <- function(d_ftb, tau) {
  curves <- fit_raw_quantile_curves(d_ftb, tau)
  hit <- make_hit_lambda_long(d_ftb, tau, curves)
  work <- hit %>%
    group_by(Participant_ID) %>%
    group_modify(~ pair_contrib_mr(.x)) %>%
    ungroup()
  sig2 <- sum(work$N) / sum(work$D)
  tibble(
    sigma2_moment = sig2,
    sigma2_moment_positive = max(0, sig2),
    n_pairs = sum(work$npair),
    n_subjects = n_distinct(d_ftb$Participant_ID),
    n_visits = nrow(d_ftb)
  )
}

run_ftb_phase_variance <- function() {
  message("Analysis 3: multi-response FTB-only phase variance")
  rows <- list()
  subject_shift_rows <- list()
  ids_ftb <- unique(dat_ftb$Participant_ID)

  for (tau in TAUS) {
    std <- standardization_from_ftb(dat_ftb)
    curves_std <- fit_std_reference_curves(dat_ftb, tau, std)
    ftb_subject_shifts <- estimate_subject_phis(dat_ftb, curves_std, tau, std, min_visits = 2L)
    ftb_subject_shifts$phi_i_centered <- ftb_subject_shifts$phi_i - mean(ftb_subject_shifts$phi_i, na.rm = TRUE)
    subject_shift_rows[[length(subject_shift_rows) + 1]] <- ftb_subject_shifts
    sig2_profile <- var(ftb_subject_shifts$phi_i_centered, na.rm = TRUE)
    boot_profile <- variance_boot_from_subject_shifts(
      ftb_subject_shifts$phi_i_centered,
      B_BOOT,
      SEED0 + 3000 + round(100 * tau)
    )
    ci_profile <- as.numeric(quantile(boot_profile, c(0.025, 0.975), na.rm = TRUE))
    p_profile_le0 <- (1 + sum(boot_profile <= 1e-10)) / (1 + length(boot_profile))

    point_moment <- estimate_mr_hit_moment(dat_ftb, tau)
    boot_moment <- rep(NA_real_, B_BOOT)
    set.seed(SEED0 + 3100 + round(100 * tau))
    for (b in seq_len(B_BOOT)) {
      db <- sample_by_subject(dat_ftb, ids_ftb)
      boot_moment[b] <- tryCatch(
        estimate_mr_hit_moment(db, tau)$sigma2_moment,
        error = function(e) NA_real_
      )
      if (b %% 50 == 0 || b == B_BOOT) {
        message("  FTB moment tau=", tau, " bootstrap ", b, "/", B_BOOT)
      }
    }
    boot_moment <- boot_moment[is.finite(boot_moment)]
    ci_moment <- if (length(boot_moment) > 10) {
      as.numeric(quantile(boot_moment, c(0.025, 0.975), na.rm = TRUE))
    } else c(NA_real_, NA_real_)
    p_moment_le0 <- if (length(boot_moment) > 10) {
      (1 + sum(boot_moment <= 0)) / (1 + length(boot_moment))
    } else NA_real_

    rows[[length(rows) + 1]] <- tibble(
      tau = tau,
      sigma2_profile = sig2_profile,
      profile_ci_low = ci_profile[1],
      profile_ci_high = ci_profile[2],
      p_profile_le0_boot = p_profile_le0,
      n_subjects_profile = nrow(ftb_subject_shifts),
      sigma2_moment = point_moment$sigma2_moment,
      sigma2_moment_positive = point_moment$sigma2_moment_positive,
      moment_ci_low = ci_moment[1],
      moment_ci_high = ci_moment[2],
      p_moment_le0_boot = p_moment_le0,
      B_moment_ok = length(boot_moment),
      n_pairs_moment = point_moment$n_pairs,
      note = "Profile variance uses FTB subject-specific composite shifts; moment variance uses cross-time multibiomarker hit covariances."
    )
  }

  res <- bind_rows(rows)
  shifts <- bind_rows(subject_shift_rows)
  write_csv(res, file.path(OUT_DIR, "analysis3_ftb_only_phase_variance.csv"))
  write_csv(shifts, file.path(OUT_DIR, "analysis3_ftb_only_profile_subject_shifts.csv"))
  res
}

# ---------------------------------------------------------------------------
# Run analyses
# ---------------------------------------------------------------------------

analysis1 <- run_mr_ga(dat)
analysis2 <- run_forward_phase()
analysis3 <- run_ftb_phase_variance()

# ---------------------------------------------------------------------------
# TeX output
# ---------------------------------------------------------------------------

make_ci <- function(lo, hi, digits = 3) {
  ifelse(is.na(lo) | is.na(hi), "--", paste0("(", fmt_num(lo, digits), ", ", fmt_num(hi, digits), ")"))
}

data_rows <- paste0(data_summary$quantity, " & ", data_summary$value, " \\\\")

ga_rows <- analysis1 %>%
  mutate(
    row = paste0(
      tau, " & ", n_biomarkers, " & ", group_df, " & ",
      fmt_num(W, 2), " & ", fmt_p(p_boot), " & ", B_ok, " \\\\"
    )
  ) %>%
  pull(row)

phase_rows <- analysis2 %>%
  mutate(
    row = paste0(
      fmt_num(tau, 2), " & ",
      fmt_num(phi_hat_weeks), " & ",
      make_ci(phi_ci_low, phi_ci_high), " & ",
      fmt_num(phi_hat_days), " & ",
      fmt_p(p_shift_boot), " & ",
      fmt_num(sigma2_phi_hat), " & ",
      make_ci(sigma2_phi_ci_low, sigma2_phi_ci_high), " & ",
      fmt_p(p_sigma2_le0_boot), " \\\\"
    )
  ) %>%
  pull(row)

ftb_rows <- analysis3 %>%
  mutate(
    row = paste0(
      fmt_num(tau, 2), " & ",
      fmt_num(sigma2_profile), " & ",
      make_ci(profile_ci_low, profile_ci_high), " & ",
      fmt_p(p_profile_le0_boot), " & ",
      fmt_num(sigma2_moment), " & ",
      fmt_num(sigma2_moment_positive), " & ",
      make_ci(moment_ci_low, moment_ci_high), " & ",
      fmt_p(p_moment_le0_boot), " & ",
      n_pairs_moment, " \\\\"
    )
  ) %>%
  pull(row)

best_ga <- analysis1 %>% arrange(p_boot) %>% slice(1)
best_phase <- analysis2 %>% arrange(p_shift_boot) %>% slice(1)

lines <- c(
  "\\section{Empirical Results for the Three Multi-response Analyses}",
  "",
  "\\subsection{Analysis Sample}",
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Sample used in the multi-response analyses. PTB is coded by \\texttt{Outcome}=1.}",
  "\\label{tab:mr_three_sample}",
  "\\small",
  "\\begin{tabular}{lr}",
  "\\toprule",
  "Quantity & Count \\\\",
  "\\midrule",
  data_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  "",
  "\\subsection{Analysis 1: Multi-response GA Spline Test}",
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Multi-response gestational-age quantile-spline tests. The tested block is the stacked PTB--FTB level and group-by-GA spline interaction across biomarkers.}",
  "\\label{tab:mr_ga_results}",
  "\\small",
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  "Quantile & Biomarkers & Group df & Wald & \\(p_{\\rm boot}\\) & Boot OK \\\\",
  "\\midrule",
  ga_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  "",
  "\\subsection{Analysis 2: Multi-response Forward-time Phase Shift}",
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Common forward-time phase-shift results. Positive \\(\\phi_\\tau\\) means PTB reaches the same FTB reference biomarker pattern earlier in gestational age. The variance column is a descriptive PTB subject-specific profile-shift variance.}",
  "\\label{tab:mr_phase_results}",
  "\\scriptsize",
  "\\begin{tabular}{lccccccc}",
  "\\toprule",
  "\\(\\tau\\) & \\(\\widehat\\phi_\\tau\\) weeks & 95\\% CI & Days & \\(p_{\\rm shift}\\) & \\(\\widehat\\sigma^2_{\\phi}\\) & 95\\% CI & \\(p_{\\le0}^{\\rm boot}\\) \\\\",
  "\\midrule",
  phase_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  "",
  "\\subsection{Analysis 3: Multi-response FTB-only Phase Variance}",
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{FTB-only multi-response phase-variance estimates. The profile estimator uses subject-specific composite phase shifts; the moment estimator uses cross-time multibiomarker quantile-hit covariances.}",
  "\\label{tab:mr_ftb_phase_variance_results}",
  "\\scriptsize",
  "\\begin{tabular}{lcccccccc}",
  "\\toprule",
  "\\(\\tau\\) & Profile \\(\\widehat\\sigma_b^2\\) & 95\\% CI & \\(p_{\\le0}\\) & Moment \\(\\widehat\\sigma_b^2\\) & Moment \\(\\widehat\\sigma_{b,+}^2\\) & 95\\% CI & \\(p_{\\le0}\\) & Pairs \\\\",
  "\\midrule",
  ftb_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  "",
  "\\subsection{Compact Readout}",
  paste0(
    "The smallest GA spline bootstrap \\(p\\)-value is ",
    fmt_p(best_ga$p_boot),
    " for the ",
    tex_escape(best_ga$tau),
    " block. "
  ),
  paste0(
    "The strongest common phase-shift evidence is at \\(\\tau=",
    fmt_num(best_phase$tau, 2),
    "\\), with \\(\\widehat\\phi_\\tau=",
    fmt_num(best_phase$phi_hat_weeks),
    "\\) weeks, or ",
    fmt_num(best_phase$phi_hat_days),
    " days. "
  ),
  "The FTB-only profile and hit-covariance analyses should be interpreted as estimates of normal timing variation, not as PTB classification analyses."
)

writeLines(lines, TEX_OUT)
message("Wrote ", TEX_OUT)
message("Wrote outputs to ", OUT_DIR)
