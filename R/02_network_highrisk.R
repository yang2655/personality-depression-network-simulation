## =========================================================
## 02_network_highrisk.R
## High-Risk profile PHQ-9 network analysis
## =========================================================

source(file.path("R", "00_setup.R"))

lpa_results <- readRDS(file.path(RDS_DIR, "lpa_results.rds"))
dat <- lpa_results$dat
full_dat4 <- lpa_results$full_dat4

dat_phq <- dat %>% dplyr::select(id, dplyr::all_of(PHQ_VARS))
full_dat_phq <- full_dat4 %>% dplyr::left_join(dat_phq, by = "id")

dat_class3_phq <- full_dat_phq %>% dplyr::filter(as.character(Class) == HIGH_RISK_CLASS)
dat_net <- dat_class3_phq %>% dplyr::select(dplyr::all_of(PHQ_VARS))

if (anyNA(dat_net)) stop("Missing PHQ-9 item data detected in the High-Risk profile.")

## ---- Main EBICglasso network ----
set.seed(SEED)
net_class3 <- bootnet::estimateNetwork(
  dat_net,
  default = "EBICglasso",
  corMethod = "cor_auto",
  tuning = 0.5
)

## ---- Figure 5: network plot with node legend ----
png(file.path(FIG_DIR, "Figure5_PHQ9_network_highrisk.png"), width = 2400, height = 1800, res = 300)
qgraph::qgraph(
  net_class3$graph,
  layout = "spring",
  labels = PHQ_VARS,
  nodeNames = unname(PHQ_LABELS),
  legend = TRUE,
  legend.cex = 0.45,
  vsize = 6,
  esize = 20,
  theme = "colorblind",
  palette = "colorblind"
)
dev.off()

## ---- Centrality and Table S3 ----
centrality_class3 <- qgraph::centrality(W)

W <- net_class3$graph
rownames(W) <- PHQ_VARS
colnames(W) <- PHQ_VARS
diag(W) <- 0

ei_values <- rowSums(W)

table_s3_ei <- data.frame(
  Node = names(ei_values),
  PHQ9_symptom = PHQ_LABELS[names(ei_values)],
  Expected_Influence = round(as.numeric(ei_values), 3),
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(dplyr::desc(Expected_Influence)) %>%
  dplyr::mutate(EI_rank = dplyr::row_number())

save_csv(table_s3_ei, "Table_S3_Expected_Influence_rankings_highrisk.csv")

## ---- Figure 4: centrality plot ----
p_centrality <- bootnet::centralityPlot(
  net_class3,
  include = c("Strength", "ExpectedInfluence"),
  orderBy = "Strength"
) +
  ggplot2::labs(x = "Centrality value", y = "PHQ-9 symptom") +
  ggplot2::scale_y_discrete(labels = PHQ_AXIS_LABELS)

ggplot2::ggsave(file.path(FIG_DIR, "Figure4_centrality_highrisk.png"),
                p_centrality, width = 8, height = 5, dpi = 600)

## ---- Bootstrap centrality stability ----
set.seed(SEED)
boot_class3 <- bootnet::bootnet(
  net_class3,
  nBoots = 2000,
  type = "case",
  statistics = c("strength", "expectedInfluence")
)

cs_values <- bootnet::corStability(boot_class3)
print(cs_values)

png(file.path(FIG_DIR, "FigureS1_case_dropping_bootstrap.png"), width = 2400, height = 1800, res = 300)
plot(boot_class3, statistics = c("strength", "expectedInfluence"), labels = TRUE, legend = TRUE)
dev.off()

set.seed(SEED)
boot_edge_class3 <- bootnet::bootnet(net_class3, nBoots = 2000, type = "nonparametric")

png(file.path(FIG_DIR, "FigureS2_edge_weight_bootstrap_CI.png"), width = 2400, height = 1800, res = 300)
plot(boot_edge_class3, labels = TRUE, order = "sample")
dev.off()

## ---- EBIC gamma sensitivity analysis ----
cor_mat <- qgraph::cor_auto(dat_net)
n_net <- nrow(dat_net)
p_net <- ncol(dat_net)
gamma_values <- c(0.25, 0.50, 1.00)

run_ebic_sensitivity <- function(gamma_value) {
  ebic_res <- qgraph::EBICglasso(S = cor_mat, n = n_net, gamma = gamma_value, returnAllResults = TRUE)

  if ("optnet" %in% names(ebic_res)) {
    W <- ebic_res$optnet
  } else if ("graph" %in% names(ebic_res)) {
    W <- ebic_res$graph
  } else {
    stop("Cannot find estimated network matrix in EBICglasso output.")
  }

  lambda_selected <- NA_real_
  if ("lambda" %in% names(ebic_res)) {
    if (length(ebic_res$lambda) == 1) {
      lambda_selected <- ebic_res$lambda
    } else if ("ebic" %in% names(ebic_res)) {
      lambda_selected <- ebic_res$lambda[which.min(ebic_res$ebic)]
    }
  }

  edge_vec <- W[upper.tri(W)]
  nonzero_edges <- edge_vec[edge_vec != 0]

  centrality_df <- data.frame(
    Node = PHQ_VARS,
    Strength = rowSums(abs(W)),
    ExpectedInfluence = rowSums(W)
  ) %>% dplyr::arrange(dplyr::desc(ExpectedInfluence))

  list(
    gamma = gamma_value,
    lambda = lambda_selected,
    W = W,
    summary = data.frame(
      gamma = gamma_value,
      lambda = lambda_selected,
      nonzero_edges = length(nonzero_edges),
      possible_edges = choose(p_net, 2),
      density = length(nonzero_edges) / choose(p_net, 2),
      global_strength = sum(abs(nonzero_edges)),
      mean_abs_edge = mean(abs(nonzero_edges)),
      min_edge = min(nonzero_edges),
      max_edge = max(nonzero_edges),
      top_EI_nodes = paste(centrality_df$Node[1:3], collapse = ", ")
    ),
    centrality = centrality_df
  )
}

sensitivity_results <- lapply(gamma_values, run_ebic_sensitivity)

gamma_sensitivity_table <- dplyr::bind_rows(lapply(sensitivity_results, function(x) x$summary))
save_csv(gamma_sensitivity_table, "Gamma_sensitivity_summary.csv")

centrality_by_gamma <- dplyr::bind_rows(
  lapply(sensitivity_results, function(x) {
    x$centrality %>%
      dplyr::mutate(gamma = x$gamma) %>%
      dplyr::group_by(gamma) %>%
      dplyr::mutate(EI_rank = rank(-ExpectedInfluence, ties.method = "min")) %>%
      dplyr::ungroup()
  })
)
save_csv(centrality_by_gamma, "Gamma_sensitivity_centrality_rankings.csv")

W_ref <- sensitivity_results[[which(gamma_values == 0.50)]]$W
edge_ref <- W_ref[upper.tri(W_ref)]
edge_ref_binary <- edge_ref != 0

edge_overlap_table <- lapply(sensitivity_results, function(x) {
  edge_temp <- x$W[upper.tri(x$W)]
  edge_temp_binary <- edge_temp != 0
  data.frame(
    gamma = x$gamma,
    edge_weight_correlation_with_gamma_0.50 = cor(edge_ref, edge_temp),
    edge_jaccard_overlap_with_gamma_0.50 =
      sum(edge_ref_binary & edge_temp_binary) / sum(edge_ref_binary | edge_temp_binary)
  )
}) %>% dplyr::bind_rows()
save_csv(edge_overlap_table, "Gamma_sensitivity_edge_overlap.csv")

## ---- Network summary statistics ----
graph_mat <- net_class3$graph
g <- igraph::graph_from_adjacency_matrix(graph_mat, mode = "undirected", weighted = TRUE, diag = FALSE)

edge_vec <- graph_mat[upper.tri(graph_mat)]
nonzero_edges <- edge_vec[edge_vec != 0]
edge_dist <- 1 / abs(igraph::E(g)$weight)

avg_path <- tryCatch(igraph::mean_distance(g, weights = edge_dist), error = function(e) igraph::mean_distance(g))
diameter_val <- tryCatch(igraph::diameter(g, weights = edge_dist), error = function(e) igraph::diameter(g))

network_summary <- data.frame(
  nodes = ncol(graph_mat),
  nonzero_edges = length(nonzero_edges),
  possible_edges = choose(ncol(graph_mat), 2),
  density = length(nonzero_edges) / choose(ncol(graph_mat), 2),
  positive_edges = sum(nonzero_edges > 0),
  negative_edges = sum(nonzero_edges < 0),
  mean_abs_edge = mean(abs(nonzero_edges)),
  min_nonzero_edge = min(nonzero_edges),
  max_nonzero_edge = max(nonzero_edges),
  average_shortest_path = avg_path,
  diameter = diameter_val,
  clustering_coefficient = igraph::transitivity(g, type = "average")
)
save_csv(network_summary, "HighRisk_network_summary.csv")

saveRDS(
  list(
    net_class3 = net_class3,
    dat_net = dat_net,
    dat_class3_phq = dat_class3_phq,
    centrality_class3 = centrality_class3,
    table_s3_ei = table_s3_ei,
    gamma_sensitivity_table = gamma_sensitivity_table,
    network_summary = network_summary
  ),
  file = file.path(RDS_DIR, "network_highrisk_results.rds")
)
