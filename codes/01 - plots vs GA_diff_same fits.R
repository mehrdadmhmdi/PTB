## All_data_spline_plots --- 

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

# keep ALL subjects (needed to plot PTB trajectories)
df <- data
#df <- df[!df$Visit_ID %in% c(5, 6, 7), ]
df <- df %>%
  mutate(
    GA = GA * 39,
    GA_delivery = GA_delivery * 39
  )

df$Outcome <- as.factor(df$Outcome)

# Time from delivery (weeks)
df$GA_diff <- df$GA - df$GA_delivery

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_qr_spline <- function(
    data, var1, var2, id,
    fit_filter = NULL,                # fit quantile curves on this subset (FTB)
    ptb_filter = NULL,                # plot PTB trajectories on this subset
    tau        = c(0.2, 0.5, 0.8),
    df_spline  = 2,
    n_grid     = 300,
    out_dir    = NULL,
    file_FTB   = "SWS_vs_GA_diff_FTB.pdf",
    file_ptb   = "SWS_vs_GA_diff_PTB.pdf",
    w = 10, h = 6
) {
  v1   <- enquo(var1); v2 <- enquo(var2)
  idq  <- enquo(id)
  fitq <- enquo(fit_filter)
  ptbq <- enquo(ptb_filter)
  
  df_all <- as_tibble(data) |>
    mutate(.x = eval_tidy(v1, data),
           .y = eval_tidy(v2, data),
           .id = eval_tidy(idq, data)) |>
    arrange(.id, .x)
  
  df_fit <- if (quo_is_null(fitq)) df_all else dplyr::filter(df_all, !!fitq)
  
  xmin <- min(df_fit$.x, na.rm = TRUE)
  xmax <- max(df_fit$.x, na.rm = TRUE)
  
  mod  <- quantreg::rq(.y ~ splines::ns(.x, df = df_spline), tau = tau, data = df_fit)
  grid <- tibble(.x = seq(xmin, xmax, length.out = n_grid))
  pred <- as.data.frame(predict(mod, newdata = grid))
  
  preds <- bind_cols(grid, pred) |>
    pivot_longer(-.x, names_to = "Quantile", values_to = "fit") |>
    mutate(Quantile = factor(gsub("^tau\\s*=\\s*", "", Quantile),
                             levels = as.character(tau)))
  
  p_FTB <- ggplot(df_fit, aes(.x, .y)) +
    geom_line(aes(group = .id), color = "grey20", alpha = .25, linewidth = .4) +
    geom_point(color = "grey20", alpha = .6, size = 1) +
    geom_line(data = preds, aes(.x, fit, color = Quantile), linewidth = 1.4) +
    labs(x = "Time From Delivery (Weeks)", y = as_label(v2), color = "Quantile") +
    scale_color_manual(
      values = c("0.2" = "#1b9e77", "0.5" = "#FF5F05", "0.8" = "#7570b3"),
      name = "Quantile")+
    theme_bw()+
    theme(legend.position = "top")
  
  df_ptb <- if (quo_is_null(ptbq)) df_all else dplyr::filter(df_all, !!ptbq)
  p_ptb <- ggplot(df_ptb, aes(.x, .y)) +
    geom_line(aes(group = .id), color = "grey30", alpha = .5, linewidth = 0.5) +
    geom_point(color = "grey30", alpha = .8, size = 1.0) +
    geom_line(data = preds, aes(.x, fit, color = Quantile), linewidth = 1.4) +
    labs(x = "Time From Delivery (Weeks)", y = as_label(v2), color = "Quantile")+
    scale_color_manual(
      values = c("0.2" = "#1b9e77", "0.5" = "#FF5F05", "0.8" = "#7570b3"),
      name = "Quantile")+
    theme_bw()+
    theme(legend.position = "top")
  
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(out_dir, file_FTB), p_FTB, width = w, height = h)
    ggsave(file.path(out_dir, file_ptb), p_ptb, width = w, height = h)
  }
  
  invisible(list(model = mod, preds = preds, plot_FTB = p_FTB, plot_ptb = p_ptb))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Single example (fit on FTB, plot on PTB)
# res <- plot_qr_spline(
#   data = df,
#   var1 = GA_diff, var2 = SWS, id = Participant_ID,
#   fit_filter = (Outcome != "1"),
#   ptb_filter = (Outcome == "1"),
#   tau = c(0.2, 0.5, 0.8),
#   df_spline = 2,
#   out_dir = "plots"
# )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5x2 grid: each ROW is one variable, COL1=FTB, COL2=PTB
vars2 <- c("SWS", "Midband", "Intercept", "AC", "Slope")

pairs_by_var <- lapply(vars2, function(v) {
  r <- plot_qr_spline(
    data = df,
    var1 = GA_diff,
    var2 = !!sym(v),
    id   = Participant_ID,
    fit_filter = (Outcome != "1"),
    ptb_filter = (Outcome == "1"),
    tau = c(0.2, 0.5, 0.8),
    df_spline = 3,
    out_dir = NULL
  )
  # one y-range per ROW (variable), shared by FTB/PTB
  ylim <- range(c(df[[v]], r$preds$fit), na.rm = TRUE)
  
  # small padding (optional but usually nicer)
  pad <- 0.05 * diff(ylim)
  if (is.finite(pad) && pad > 0) ylim <- ylim + c(-pad, pad)
  
  left  <- r$plot_FTB + ggtitle("FTB Group")  + coord_cartesian(ylim = ylim) +
    theme(legend.position = "bottom")
  right <- r$plot_ptb + ggtitle("PTB Group")  + coord_cartesian(ylim = ylim) +
    theme(legend.position = "bottom")
  
  left | right
})

grid_5x2 <- wrap_plots(pairs_by_var, ncol = 1)
grid_5x2
dir.create("plots", showWarnings = FALSE, recursive = TRUE)
ggsave("plots/GA_diff_grid_5x2.pdf", grid_5x2, width = 12, height = 24)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2x5 grid (row1=FTB, row2=PTB), optional
plots_FTB <- lapply(vars2, function(v) {
  plot_qr_spline(df, GA_diff, !!sym(v), Participant_ID,
                 fit_filter = (Outcome != "1"),
                 ptb_filter = (Outcome == "1"),
                 df_spline = 3, out_dir = NULL)$plot_FTB +
    ggtitle("FTB Group") + theme(legend.position = "top")
})
plots_ptb <- lapply(vars2, function(v) {
  plot_qr_spline(df, GA_diff, !!sym(v), Participant_ID,
                 fit_filter = (Outcome != "1"),
                 ptb_filter = (Outcome == "1"),
                 df_spline = 3, out_dir = NULL)$plot_ptb +
    ggtitle("PTB Group")  + theme(legend.position = "top")
})

grid_2x5 <- wrap_plots(plots_FTB, ncol = 5) / wrap_plots(plots_ptb, ncol = 5)
grid_2x5
ggsave("plots/GA_diff_grid_2x5.pdf", grid_2x5, width = 24, height = 8)
