## =========================================================
## 04_insilico_simulation.R
## Exploratory in silico simulation-based perturbation analysis
## =========================================================

source(file.path("R", "00_setup.R"))

network_results <- readRDS(file.path(RDS_DIR, "network_highrisk_results.rds"))
net_class3 <- network_results$net_class3
dat_net <- network_results$dat_net

graph_mat <- if (is.list(net_class3) && !is.null(net_class3$graph)) net_class3$graph else net_class3
stopifnot(is.matrix(graph_mat))

## ---- Expected Influence order ----
ei <- rowSums(graph_mat)
names(ei) <- colnames(graph_mat)
order_nodes <- names(sort(ei, decreasing = TRUE))

## ---- Gaussian conditional simulation parameters ----
mu <- colMeans(dat_net, na.rm = TRUE)
Sigma <- stats::cov(dat_net, use = "pairwise.complete.obs")
TL_baseline <- sum(mu)

do_TL_gaussian <- function(targets, k = -1, mu, Sigma) {
  all_vars <- names(mu)
  if (!all(targets %in% all_vars)) {
    stop("Some target variables are not present in mu: ", paste(setdiff(targets, all_vars), collapse = ", "))
  }

  idx_t <- match(targets, all_vars)
  idx_r <- setdiff(seq_along(all_vars), idx_t)

  mu_t <- mu[idx_t]
  mu_r <- mu[idx_r]
  Sigma_tt <- Sigma[idx_t, idx_t, drop = FALSE]
  Sigma_rt <- Sigma[idx_r, idx_t, drop = FALSE]

  if (any(is.na(Sigma_tt)) || det(Sigma_tt) == 0) {
    stop("Sigma_tt is singular or contains missing values.")
  }

  ## Target nodes after perturbation: mean shifted by k SD.
  ## Values are truncated to the theoretical PHQ-9 item range of 0-3.
  sd_t <- sqrt(diag(Sigma_tt))
  x_star <- mu_t + k * sd_t
  x_star <- pmin(3, pmax(0, x_star))

  if (length(idx_r) > 0) {
    mu_r_given <- as.numeric(mu_r + Sigma_rt %*% solve(Sigma_tt) %*% (x_star - mu_t))
    mu_r_given <- pmin(3, pmax(0, mu_r_given))
    TL_after <- sum(x_star) + sum(mu_r_given)
  } else {
    TL_after <- sum(x_star)
  }

  list(
    TL_after = TL_after,
    dTL = TL_after - sum(mu),
    dTL_pct = (TL_after - sum(mu)) / sum(mu) * 100
  )
}

k_shift <- -1
n_nodes <- length(order_nodes)

## ---- Cumulative perturbations in EI order ----
delta_TL <- numeric(n_nodes)
delta_TL_pct <- numeric(n_nodes)

for (m in seq_len(n_nodes)) {
  res_m <- do_TL_gaussian(targets = order_nodes[1:m], k = k_shift, mu = mu, Sigma = Sigma)
  delta_TL[m] <- res_m$dTL
  delta_TL_pct[m] <- res_m$dTL_pct
}

## ---- Single-node perturbations ----
single_delta_TL <- numeric(n_nodes)
single_delta_TL_pct <- numeric(n_nodes)

for (i in seq_along(order_nodes)) {
  res_i <- do_TL_gaussian(targets = order_nodes[i], k = k_shift, mu = mu, Sigma = Sigma)
  single_delta_TL[i] <- res_i$dTL
  single_delta_TL_pct[i] <- res_i$dTL_pct
}

simulation_results <- data.frame(
  node = order_nodes,
  symptom = PHQ_LABELS[order_nodes],
  EI = as.numeric(ei[order_nodes]),
  single_delta_TL = single_delta_TL,
  single_delta_TL_pct = single_delta_TL_pct,
  cumulative_delta_TL = delta_TL,
  cumulative_delta_TL_pct = delta_TL_pct,
  stringsAsFactors = FALSE
)

print(simulation_results)
save_csv(simulation_results, "Simulation_results_highrisk.csv")

## ---- Figure 6 ----
df_plot <- data.frame(
  node = factor(order_nodes, levels = order_nodes),
  single_pct = single_delta_TL_pct,
  cum_pct = delta_TL_pct
)

y_min <- min(df_plot$single_pct, df_plot$cum_pct) - 5

p_sim <- ggplot2::ggplot(df_plot, ggplot2::aes(x = node)) +
  ggplot2::geom_col(ggplot2::aes(y = single_pct, fill = "Single-node"), width = 0.6) +
  ggplot2::geom_text(ggplot2::aes(y = single_pct, label = sprintf("%.1f%%", single_pct)),
                     vjust = -0.5, size = 3.8) +
  ggplot2::geom_line(ggplot2::aes(y = cum_pct, color = "Cumulative", group = 1), linewidth = 1) +
  ggplot2::geom_point(ggplot2::aes(y = cum_pct, color = "Cumulative"), size = 2.3) +
  ggplot2::geom_text(data = df_plot[-1, ],
                     ggplot2::aes(y = cum_pct, label = sprintf("%.1f%%", cum_pct)),
                     vjust = -1.2, size = 3.5, color = "black") +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::scale_fill_manual(values = c("Single-node" = "grey70")) +
  ggplot2::scale_color_manual(values = c("Cumulative" = "black")) +
  ggplot2::scale_y_continuous(
    limits = c(y_min, 5),
    breaks = seq(floor(y_min / 10) * 10, 0, by = 10),
    labels = function(x) paste0(x, "%")
  ) +
  ggplot2::labs(x = "PHQ-9 symptoms (EI order)", y = expression(Delta * TL ~ "(%)")) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = ggplot2::element_blank()
  )

ggplot2::ggsave(file.path(FIG_DIR, "Figure6_in_silico_simulation.png"),
                p_sim, width = 8, height = 5, dpi = 600)

saveRDS(
  list(simulation_results = simulation_results, TL_baseline = TL_baseline,
       k_shift = k_shift, order_nodes = order_nodes),
  file = file.path(RDS_DIR, "simulation_results.rds")
)
