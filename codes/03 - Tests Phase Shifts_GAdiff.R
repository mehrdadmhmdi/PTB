#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Phase-shift (horizontal shift) quantile regression on delivery-aligned time x = GA_diff
# Model: Q_tau(Y | x, G) = alpha_tau + delta_tau * G + f_tau(x + phi_tau * G),
# where G=1 for PTB, 0 for FTB; f_tau is a natural spline.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
library(readr)
library(dplyr)
library(tidyr)
library(quantreg)
library(splines)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- read_csv("../../data/V1_V7_longformat_k1.csv")

df <- data
df <- df %>% mutate(GA = GA * 39, GA_delivery = GA_delivery * 39)
df$Outcome <- as.factor(df$Outcome)
df$GA_diff <- df$GA - df$GA_delivery
# -----------------------------
# Helpers
# -----------------------------
qr_checkloss <- function(mod, tau) {
  u <- resid(mod)
  sum(u * (tau - (u < 0)))
}

# Fit rq for a given phi and return objective (check-loss)
fit_at_phi <- function(dat, yvar, tau, df_spline, phi, G_col = "G", x_col = "GA_diff") {
  dat <- dat %>% mutate(x_shift = .data[[x_col]] + phi * .data[[G_col]])
  f <- as.formula(sprintf("%s ~ ns(x_shift, df=%d)", yvar, df_spline))
  m <- tryCatch(rq(f, tau = tau, data = dat), error = function(e) NULL)
  if (is.null(m)) return(list(obj = NA_real_, mod = NULL))
  list(obj = qr_checkloss(m, tau), mod = m)
}

# Profile over a grid of phi and return the minimiser + profile table
profile_phi <- function(dat, yvar, tau, df_spline = 3,
                        phi_grid = seq(-6, 6, by = 0.1),
                        G_col = "G", x_col = "GA_diff") {
  
  objs <- rep(NA_real_, length(phi_grid))
  for (k in seq_along(phi_grid)) {
    objs[k] <- fit_at_phi(dat, yvar, tau, df_spline, phi_grid[k], G_col, x_col)$obj
  }
  
  ok <- is.finite(objs)
  if (!any(ok)) {
    return(list(phi_hat = NA_real_, obj_min = NA_real_, obj_phi0 = NA_real_,
                prof = data.frame(phi = phi_grid, obj = objs)))
  }
  
  objs2 <- objs
  objs2[!ok] <- Inf          # force NA/Inf to never be chosen
  k_min <- which.min(objs2)
  phi_hat <- phi_grid[k_min]
  obj_min <- objs2[k_min]
  
  # objective at phi=0 (interpolate if 0 not in grid)
  if (any(phi_grid == 0)) {
    obj_phi0 <- objs[which(phi_grid == 0)[1]]
  } else {
    # linear interpolation as fallback
    obj_phi0 <- approx(phi_grid[ok], objs[ok], xout = 0, rule = 2)$y
  }
  
  list(
    phi_hat = phi_hat,
    obj_min = obj_min,
    obj_phi0 = obj_phi0,
    prof = data.frame(phi = phi_grid, obj = objs)
  )
}




# ============================================================
# Random-effect phase shift: phi_i = phi + b_i,  b_i ~ N(0, sigma^2)
# Model used in fitting: Q_tau(Y|x,G,ID)= alpha + delta*G + f(x + (phi + b_ID)*G)
# We estimate (phi, {b_ID}) by alternating updates, and report sigma^2 = var(b_ID).
# ============================================================

rho_tau_sum <- function(u, tau) sum(u * (tau - (u < 0)))

# Fit rq given current phi and b_by_id (named vector for PTB subjects; 0 for FTB)
fit_given_phi_b <- function(dat, yvar, tau, df_spline, phi, b_by_id,
                            id_col="Participant_ID", G_col="G", x_col="GA_diff") {
  dat <- dat
  bid <- as.character(dat[[id_col]])
  bvec <- ifelse(dat[[G_col]] == 1L, b_by_id[bid], 0)
  bvec[is.na(bvec)] <- 0
  dat$x_shift <- dat[[x_col]] + (phi + bvec) * dat[[G_col]]
  
  f <- as.formula(sprintf("%s ~  ns(x_shift, df=%d)", yvar, df_spline))
  m <- tryCatch(rq(f, tau = tau, data = dat), error = function(e) NULL)
  if (is.null(m)) return(list(obj = Inf, mod = NULL, dat = dat))
  list(obj = qr_checkloss(m, tau), mod = m, dat = dat)
}

# Alternating minimization for fixed sigma (shrink b_i via N(0, sigma^2))
fit_phase_shift_RE_given_sigma <- function(dat, yvar, tau,
                                           id_col="Participant_ID",
                                           x_col="GA_diff",
                                           G_col="G",
                                           df_spline=3,
                                           phi_grid=seq(-10,10,by=0.25),
                                           b_grid=seq(-6,6,by=0.25),
                                           sigma=1,
                                           max_iter=10,
                                           tol=1e-3,
                                           seed=1) {
  set.seed(seed)
  
  # PTB subject IDs
  G_by_id <- tapply(dat[[G_col]], dat[[id_col]], function(z) z[1])
  ptb_ids <- names(G_by_id)[G_by_id == 1L]
  
  # initialize
  b_by_id <- setNames(rep(0, length(ptb_ids)), ptb_ids)
  phi <- 0
  
  # helper: enforce mean(b)=0 by absorbing into phi
  recenter_b <- function(phi, b) {
    mb <- mean(b, na.rm = TRUE)
    list(phi = phi + mb, b = b - mb)
  }
  
  fit_best <- NULL
  
  for (it in seq_len(max_iter)) {
    phi_old <- phi
    b_old   <- b_by_id
    
    # --- (1) update phi by grid search, holding b fixed (refit rq each grid point) ---
    objs_phi <- rep(Inf, length(phi_grid))
    mods_phi <- vector("list", length(phi_grid))
    for (k in seq_along(phi_grid)) {
      tmp <- fit_given_phi_b(dat, yvar, tau, df_spline, phi_grid[k], b_by_id,
                             id_col=id_col, G_col=G_col, x_col=x_col)
      objs_phi[k] <- tmp$obj
      mods_phi[[k]] <- tmp
    }
    kmin <- which.min(objs_phi)
    phi  <- phi_grid[kmin]
    fit_best <- mods_phi[[kmin]]
    mod_cur  <- fit_best$mod
    
    # --- (2) update each PTB subject b_i by 1D grid, holding mod_cur fixed ---
    for (sid in ptb_ids) {
      rows_i <- which(as.character(dat[[id_col]]) == sid & dat[[G_col]] == 1L)
      if (length(rows_i) == 0) next
      
      dat_i <- dat[rows_i, , drop = FALSE]
      y_i   <- dat_i[[yvar]]
      
      # evaluate candidate b values
      loss_cand <- rep(Inf, length(b_grid))
      for (kk in seq_along(b_grid)) {
        b_cand <- b_grid[kk]
        dat_i2 <- dat_i
        dat_i2$x_shift <- dat_i2[[x_col]] + (phi + b_cand)  # because G=1 inside PTB
        yhat <- tryCatch(predict(mod_cur, newdata = dat_i2), error = function(e) rep(NA_real_, length(y_i)))
        if (anyNA(yhat)) next
        u <- y_i - yhat
        loss_cand[kk] <- rho_tau_sum(u, tau) + 0.5 * (b_cand^2) / (sigma^2)
      }
      b_by_id[sid] <- b_grid[which.min(loss_cand)]
    }
    
    # recenter b so mean(b)=0, absorb mean into phi (identifiability)
    tmpc <- recenter_b(phi, b_by_id)
    phi <- tmpc$phi
    b_by_id <- tmpc$b
    
    # convergence check
    dphi <- abs(phi - phi_old)
    db   <- max(abs(b_by_id - b_old), na.rm = TRUE)
    if (max(dphi, db) < tol) break
  }
  
  # final refit using final (phi, b)
  fit_fin <- fit_given_phi_b(dat, yvar, tau, df_spline, phi, b_by_id,
                             id_col=id_col, G_col=G_col, x_col=x_col)
  
  # sigma^2 estimate from BLUP-like b's (random-effect scale summary)
  sig2_hat <- var(as.numeric(b_by_id), na.rm = TRUE)
  
  list(phi_hat = phi, b_hat = b_by_id, sig2_hat = sig2_hat,
       obj = fit_fin$obj, mod = fit_fin$mod)
}

# choose sigma by grid-search (simple, stable); includes log(sigma) complexity term
fit_phase_shift_RE <- function(dat, yvar, tau,
                               id_col="Participant_ID",
                               x_col="GA_diff",
                               G_col="G",
                               df_spline=3,
                               phi_grid=seq(-10,10,by=0.25),
                               b_grid=seq(-6,6,by=0.25),
                               sigma_grid=c(0.25,0.5,1,2,4),
                               max_iter=10,
                               tol=1e-3,
                               seed=1) {
  
  # PTB count for the complexity term
  G_by_id <- tapply(dat[[G_col]], dat[[id_col]], function(z) z[1])
  n_ptb   <- sum(G_by_id == 1L, na.rm = TRUE)
  
  best <- NULL
  best_crit <- Inf
  
  for (sg in sigma_grid) {
    fit <- fit_phase_shift_RE_given_sigma(dat, yvar, tau,
                                          id_col=id_col, x_col=x_col, G_col=G_col,
                                          df_spline=df_spline, phi_grid=phi_grid, b_grid=b_grid,
                                          sigma=sg, max_iter=max_iter, tol=tol, seed=seed + round(100*sg))
    # pseudo-REML-ish criterion to prevent sigma -> infinity overfit
    crit <- fit$obj + 0.5 * sum(fit$b_hat^2, na.rm = TRUE) / (sg^2) + n_ptb * log(sg)
    
    if (crit < best_crit) {
      best_crit <- crit
      best <- fit
      best$sigma_hat <- sg
      best$crit <- crit
    }
  }
  best
}




# -----------------------------
# Cluster bootstrap for phi_hat and a simple "improvement over phi=0" statistic
# T = Obj(phi=0) - min_phi Obj(phi); larger T = more evidence for a nonzero shift
# -----------------------------
# --- map ptb_level robustly (accept "PTB" or numeric index) ---
resolve_ptb_level <- function(dat, outcome_col, ptb_level){
  levs <- levels(dat[[outcome_col]])
  if (is.numeric(ptb_level) && length(ptb_level) == 1L) {
    if (ptb_level >= 1 && ptb_level <= length(levs)) return(levs[ptb_level])
  }
  if (is.character(ptb_level) && ptb_level %in% levs) return(ptb_level)
  stop("ptb_level must be a valid Outcome level name (e.g., 'PTB') or a level index.")
}

# -----------------------------
# Cluster bootstrap CI for phi_hat + PERMUTATION test p-value for T
# -----------------------------
phase_shift_cluster_boot <- function(df, yvar, tau,
                                     id_col = "Participant_ID",
                                     outcome_col = "Outcome",
                                     ptb_level = 1,
                                     x_col = "GA_diff",
                                     G_col = "G",
                                     df_spline = 3,
                                     phi_grid = seq(-10, 10, by = 0.25),
                                     b_grid   = seq(-6, 6, by = 0.25),
                                     sigma_grid = c(0.25, 0.5, 1, 2, 4),
                                     B = 100,
                                     B_sigma0 = 100,
                                     seed = 2026) {
  
  dat <- df
  dat[[id_col]]      <- as.factor(dat[[id_col]])
  dat[[outcome_col]] <- as.factor(dat[[outcome_col]])
  
  ptb_lab <- resolve_ptb_level(dat, outcome_col, ptb_level)
  dat <- dat %>% mutate(G = as.integer(.data[[outcome_col]] == ptb_lab))
  
  # cluster index
  idx_by_id <- split(seq_len(nrow(dat)), dat[[id_col]])
  uniq_ids  <- names(idx_by_id)
  n_ids     <- length(uniq_ids)
  
  # ---------- observed (RE) ----------
  fit_re <- fit_phase_shift_RE(dat, yvar, tau,
                               id_col=id_col, x_col=x_col, G_col=G_col,
                               df_spline=df_spline,
                               phi_grid=phi_grid, b_grid=b_grid,
                               sigma_grid=sigma_grid,
                               seed=seed)
  
  phi_hat_obs  <- fit_re$phi_hat
  sig2_hat_obs <- fit_re$sig2_hat
  
  # ---------- null sigma^2 = 0 (fixed shift) ----------
  obs0 <- profile_phi(dat, yvar, tau, df_spline, phi_grid, G_col = G_col, x_col = x_col)
  obj0 <- obs0$obj_min
  
  # RE improvement statistic (bigger => more evidence sigma^2>0)
  Tsig_obs <- obj0 - fit_re$obj
  
  # ---------- cluster bootstrap CI for (phi, sigma^2) under the fitted RE procedure ----------
  set.seed(seed + 10)
  phi_boot  <- rep(NA_real_, B)
  sig2_boot <- rep(NA_real_, B)
  Tsig_boot <- rep(NA_real_, B)
  
  for (b in seq_len(B)) {
    samp_ids <- sample(uniq_ids, size = n_ids, replace = TRUE)
    rows_b   <- unlist(idx_by_id[samp_ids], use.names = FALSE)
    dat_b    <- dat[rows_b, , drop = FALSE]
    
    fb <- fit_phase_shift_RE(dat_b, yvar, tau,
                             id_col=id_col, x_col=x_col, G_col=G_col,
                             df_spline=df_spline,
                             phi_grid=phi_grid, b_grid=b_grid,
                             sigma_grid=sigma_grid,
                             seed=seed + 1000 + b)
    
    phi_boot[b]  <- fb$phi_hat
    sig2_boot[b] <- fb$sig2_hat
    
    # null-fit on bootstrap sample
    b0 <- profile_phi(dat_b, yvar, tau, df_spline, phi_grid, G_col = G_col, x_col = x_col)
    Tsig_boot[b] <- b0$obj_min - fb$obj
  }
  
  ci_phi  <- quantile(phi_boot,  probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  ci_sig2 <- quantile(sig2_boot, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  
  # ---------- test sigma^2=0 via wild (cluster) bootstrap under H0 ----------
  # fit H0 model at phi_hat0 and generate y* = yhat0 + w_id * resid0 (w_id = +/-1)
  set.seed(seed + 999)
  # build H0 fitted model at phi_hat0
  dat0 <- dat %>% mutate(x_shift = .data[[x_col]] + obs0$phi_hat * .data[[G_col]])
  f0 <- as.formula(sprintf("%s ~ ns(x_shift, df=%d)", yvar, df_spline))
  m0 <- rq(f0, tau = tau, data = dat0)
  yhat0 <- predict(m0, newdata = dat0)
  r0 <- dat0[[yvar]] - yhat0
  
  # subject weights
  w_id <- setNames(sample(c(-1, 1), size = length(uniq_ids), replace = TRUE), uniq_ids)
  
  Tsig_null <- rep(NA_real_, B_sigma0)
  for (bb in seq_len(B_sigma0)) {
    w_id <- setNames(sample(c(-1, 1), size = length(uniq_ids), replace = TRUE), uniq_ids)
    dat_star <- dat0
    dat_star[[yvar]] <- yhat0 + w_id[as.character(dat0[[id_col]])] * r0
    
    # refit RE on y*
    fb <- fit_phase_shift_RE(dat_star, yvar, tau,
                             id_col=id_col, x_col=x_col, G_col=G_col,
                             df_spline=df_spline,
                             phi_grid=phi_grid, b_grid=b_grid,
                             sigma_grid=sigma_grid,
                             seed=seed + 5000 + bb)
    
    # refit H0 on y*
    b0 <- profile_phi(dat_star, yvar, tau, df_spline, phi_grid, G_col = G_col, x_col = x_col)
    Tsig_null[bb] <- b0$obj_min - fb$obj
  }
  
  p_sigma0 <- mean(Tsig_null >= Tsig_obs, na.rm = TRUE)
  
  list(
    yvar = yvar, tau = tau,
    phi_hat = phi_hat_obs,
    phi_ci_low = ci_phi[1], phi_ci_high = ci_phi[2],
    sig2_hat = sig2_hat_obs,
    sig2_ci_low = ci_sig2[1], sig2_ci_high = ci_sig2[2],
    Tsig_obs = Tsig_obs,
    p_sigma0 = p_sigma0
  )
}

# -----------------------------
# Example usage
# -----------------------------
# # Example for one biomarker and tau:
 # res <- phase_shift_cluster_boot(df, yvar="SWS", tau=0.5,
 #                                id_col="Participant_ID",
 #                                outcome_col="Outcome",
 #                                ptb_level=1,
 #                                x_col="GA_diff",
 #                                G_col ="G",
 #                                df_spline=3,
 #                                phi_grid=seq(-10,10,by=0.1),
 #                                B=100, seed=2026)
 # 
 # print(res[c("yvar","tau","phi_hat","phi_ci_low","phi_ci_high","p_perm_T")])
 # ## Profile plot (objective vs phi)
 # plot(res$prof$phi, res$prof$obj, type="l", xlab="phi (weeks)", ylab="check-loss")
 # abline(v=0, lty=2)
 # abline(v=res$phi_hat, lty=3)

# -----------------------------
# Run for multiple biomarkers/taus and export a summary table
# -----------------------------
run_phase_shift_grid <- function(df, vars, taus,
                                 id_col="Participant_ID",
                                 outcome_col="Outcome",
                                 ptb_level=1,
                                 x_col="GA_diff",
                                 G_col ="G",
                                 df_spline=3,
                                 phi_grid=seq(-10,10,by=0.5),
                                 B=100, seed=2026) {
  
  out <- list()
  k <- 0
  for (v in vars) {
    for (t in taus) {
      k <- k + 1
      res <- phase_shift_cluster_boot(df, yvar=v, tau=t,
                                      id_col=id_col,
                                      outcome_col=outcome_col,
                                      ptb_level=ptb_level,
                                      x_col=x_col,
                                      G_col = G_col,
                                      df_spline=df_spline,
                                      phi_grid=phi_grid,
                                      B=B, seed=seed + k)
      out[[k]] <- data.frame(
        biomarker = v, tau = t,
        phi_hat = res$phi_hat,
        phi_ci_low = res$phi_ci_low,
        phi_ci_high = res$phi_ci_high,
        sig2_hat = res$sig2_hat,
        sig2_ci_low = res$sig2_ci_low,
        sig2_ci_high = res$sig2_ci_high,
        Tsig_obs = res$Tsig_obs,
        p_sigma0 = res$p_sigma0
      )
      
    }
  }
  bind_rows(out) %>% arrange(biomarker, tau)
}

# Example:
vars2 <- c("SWS","Midband","Intercept","AC","Slope")
taus  <- c(0.2,0.5,0.8)


summ <- run_phase_shift_grid(df, vars2, taus, ptb_level=1, B=100)
write.csv(summ, "Test_results/qr_phase_shift_summary_GAdiff.csv", row.names=FALSE)

