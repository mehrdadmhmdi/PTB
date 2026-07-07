## Group-wise spline QR plots: fit separately on FTB and PTB
rm(list = ls())

library(readr)
library(dplyr)
library(tidyr)
library(quantreg)
library(ggplot2)
library(splines)
library(rlang)
library(patchwork)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- read_csv("../../data/V1_V7_longformat_k1.csv")

df <- data %>%
  mutate(
    GA = GA * 39,
    GA_delivery = GA_delivery * 39,
    Outcome = as.factor(Outcome),
    GA_diff = GA - GA_delivery
  )

df_FTB <- df %>% filter(Outcome != "1")
df_PTB <- df %>% filter(Outcome == "1")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit_qr_spline <- function(data, var1, var2, id,
                          tau = c(0.2, 0.5, 0.8),
                          df_spline = 3,
                          n_grid = 300) {
  v1  <- enquo(var1); v2 <- enquo(var2); idq <- enquo(id)
  
  d <- as_tibble(data) %>%
    mutate(.x = !!v1, .y = !!v2, .id = !!idq) %>%
    filter(!is.na(.x), !is.na(.y), !is.na(.id)) %>%
    arrange(.id, .x)
  
  xmin <- min(d$.x, na.rm = TRUE)
  xmax <- max(d$.x, na.rm = TRUE)
  
  mod  <- quantreg::rq(.y ~ splines::ns(.x, df = df_spline),
                       tau = tau, data = d)
  
  grid <- tibble(.x = seq(xmin, xmax, length.out = n_grid))
  pred <- as.data.frame(predict(mod, newdata = grid))
  
  preds <- bind_cols(grid, pred) %>%
    pivot_longer(-.x, names_to = "Quantile", values_to = "fit") %>%
    mutate(
      Quantile = factor(gsub("^tau\\s*=\\s*", "", Quantile),
                        levels = as.character(tau))
    )
  
  list(data = d, model = mod, preds = preds, x_rng = c(xmin, xmax))
}

make_plot <- function(fit_obj, title, ylab,
                      xlim = NULL, ylim = NULL,
                      show_legend = TRUE) {
  p <- ggplot(fit_obj$data, aes(.x, .y)) +
    geom_line(aes(group = .id), color = "grey30", alpha = 0.35, linewidth = 0.4) +
    geom_point(color = "grey30", alpha = 0.70, size = 0.9) +
    geom_line(data = fit_obj$preds, aes(.x, fit, color = Quantile), linewidth = 1.2) +
    labs(x = "Time From Delivery (Weeks)", y = ylab, color = "Quantile", title = title) +
    scale_color_manual(
      values = c("0.2" = "#1b9e77", "0.5" = "#FF5F05", "0.8" = "#7570b3"),
      name = "Quantile"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (!show_legend) p <- p + theme(legend.position = "none")
  if (!is.null(xlim) || !is.null(ylim)) p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  p
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vars2 <- c("SWS", "Midband", "Intercept", "AC", "Slope")

# Pre-fit once per biomarker (FTB + PTB)
fits <- lapply(vars2, function(v) {
  f_ftb <- fit_qr_spline(df_FTB, GA_diff, !!sym(v), Participant_ID,
                         tau = c(0.2, 0.5, 0.8), df_spline = 3, n_grid = 300)
  f_ptb <- fit_qr_spline(df_PTB, GA_diff, !!sym(v), Participant_ID,
                         tau = c(0.2, 0.5, 0.8), df_spline = 3, n_grid = 300)
  
  # shared x/y limits for comparability (per biomarker)
  xlim_all <- range(c(f_ftb$data$.x, f_ptb$data$.x), na.rm = TRUE)
  ylim_all <- range(c(f_ftb$data$.y, f_ptb$data$.y,
                      f_ftb$preds$fit, f_ptb$preds$fit), na.rm = TRUE)
  
  pad <- 0.05 * diff(ylim_all)
  if (is.finite(pad) && pad > 0) ylim_all <- ylim_all + c(-pad, pad)
  
  list(v = v, ftb = f_ftb, ptb = f_ptb, xlim = xlim_all, ylim = ylim_all)
})
names(fits) <- vars2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5x2 grid: each ROW is one biomarker, COL1=FTB, COL2=PTB
pairs_by_var <- lapply(vars2, function(v) {
  obj <- fits[[v]]
  
  p1 <- make_plot(obj$ftb, "FTB Group", v, xlim = obj$xlim, ylim = obj$ylim, show_legend = TRUE)
  p2 <- make_plot(obj$ptb, "PTB Group", v, xlim = obj$xlim, ylim = obj$ylim, show_legend = TRUE)
  
  p1 | p2
})

grid_5x2 <- wrap_plots(pairs_by_var, ncol = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

dir.create("plots", showWarnings = FALSE, recursive = TRUE)
ggsave("plots/GA_diff_grid_5x2.pdf", grid_5x2, width = 12, height = 24)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2x5 grid: row1=FTB, row2=PTB
plots_FTB <- lapply(vars2, function(v) {
  obj <- fits[[v]]
  make_plot(obj$ftb, "FTB Group", v, xlim = obj$xlim, ylim = obj$ylim, show_legend = TRUE)
})

plots_PTB <- lapply(vars2, function(v) {
  obj <- fits[[v]]
  make_plot(obj$ptb, "PTB Group", v, xlim = obj$xlim, ylim = obj$ylim, show_legend = TRUE)
})

grid_2x5 <- (wrap_plots(plots_FTB, ncol = 5) /
               wrap_plots(plots_PTB, ncol = 5)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("plots/GA_diff_grid_2x5.pdf", grid_2x5, width = 24, height = 8)
