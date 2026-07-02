## =========================================================
## 01_lpa_profiles.R
## Latent Profile Analysis and class-level PHQ-9 comparisons
## =========================================================

source(file.path("R", "00_setup.R"))

dat <- load_clean_data()
lpa_dat <- make_lpa_data(dat)

## ---- LPA: 1- to 5-class solutions, tidyLPA Model 1 ----
set.seed(SEED)
model_list <- tidyLPA::estimate_profiles(
  lpa_dat[, c("E", "N", "P", "L")],
  n_profiles = 1:5,
  models = 1,
  seed = SEED,
  nrep = 5
)

fit_info <- tidyLPA::get_fit(model_list)
print(fit_info)
save_csv(fit_info, "Table_S1_LPA_fit_indices.csv")

## ---- Select 4-class solution ----
idx4 <- which(fit_info$Classes == 4)
model4 <- model_list[[idx4]]
dat4 <- tidyLPA::get_data(model4)
dat4$id <- lpa_dat$id

full_dat4 <- lpa_dat %>%
  dplyr::left_join(dat4[, c("id", "Class")], by = "id") %>%
  dplyr::mutate(Class = factor(Class))

print(table(full_dat4$Class))
print(prop.table(table(full_dat4$Class)))

## ---- Profile summary: full-sample z scores ----
dat4_z <- full_dat4 %>%
  dplyr::mutate(
    E_z = as.numeric(scale(E)),
    N_z = as.numeric(scale(N)),
    P_z = as.numeric(scale(P)),
    L_z = as.numeric(scale(L))
  )

dat4_long <- dat4_z %>%
  dplyr::select(Class, E_z, N_z, P_z, L_z) %>%
  tidyr::pivot_longer(
    cols = c(E_z, N_z, P_z, L_z),
    names_to = "Dimension",
    values_to = "Zscore"
  ) %>%
  dplyr::mutate(
    Dimension = factor(
      Dimension,
      levels = c("E_z", "N_z", "P_z", "L_z"),
      labels = c("E (Extraversion)", "N (Neuroticism)", "P (Psychoticism)", "L (Lie)")
    ),
    Class = factor(Class)
  )

profile_summary <- dat4_long %>%
  dplyr::group_by(Class, Dimension) %>%
  dplyr::summarise(
    mean_z = mean(Zscore, na.rm = TRUE),
    sd_z = sd(Zscore, na.rm = TRUE),
    n = dplyr::n(),
    se_z = sd_z / sqrt(n),
    .groups = "drop"
  )

save_csv(profile_summary, "Table_S2_profile_mean_z_scores.csv")

p_profile <- ggplot2::ggplot(
  profile_summary,
  ggplot2::aes(x = Dimension, y = mean_z, group = Class, color = Class)
) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
    width = .08,
    linewidth = .8
  ) +
  ggplot2::scale_color_brewer(palette = "Set1") +
  ggplot2::labs(x = "EPQ Dimensions", y = "Standardized score (z-score)", color = "Class") +
  theme_manuscript(base_size = 16)

ggplot2::ggsave(
  file.path(FIG_DIR, "Figure2_EPQ_profile_plot.png"),
  p_profile,
  width = 8,
  height = 5,
  dpi = 600
)

## ---- PHQ-9 descriptive statistics by class ----
N_total <- nrow(full_dat4)

table1_panelA <- full_dat4 %>%
  dplyr::mutate(
    Class = as.character(Class),
    PHQ9_total = DEP_total,
    PHQ9_ge10 = PHQ9_total >= 10
  ) %>%
  dplyr::group_by(Class) %>%
  dplyr::summarise(
    n = dplyr::n(),
    PHQ9_mean = mean(PHQ9_total, na.rm = TRUE),
    PHQ9_sd = sd(PHQ9_total, na.rm = TRUE),
    PHQ9_ge10_n = sum(PHQ9_ge10, na.rm = TRUE),
    PHQ9_ge10_percent = 100 * mean(PHQ9_ge10, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(percent = 100 * n / N_total) %>%
  dplyr::left_join(CLASS_LABELS, by = "Class") %>%
  dplyr::arrange(as.numeric(Class))

save_csv(table1_panelA, "Table1_PanelA_PHQ9_descriptives_by_class.csv")

fit_aov_dep <- stats::aov(DEP_total ~ Class, data = full_dat4)
print(summary(fit_aov_dep))
print(effectsize::eta_squared(fit_aov_dep, partial = FALSE))

tuk_dep <- stats::TukeyHSD(fit_aov_dep)
print(tuk_dep)

tukey_table <- as.data.frame(tuk_dep$Class)
tukey_table$Comparison <- rownames(tukey_table)
rownames(tukey_table) <- NULL
tukey_table <- tukey_table %>% dplyr::select(Comparison, diff, lwr, upr, `p adj`)
save_csv(tukey_table, "Table1_PanelB_TukeyHSD_PHQ9.csv")

p_dep_box <- ggplot2::ggplot(full_dat4, ggplot2::aes(x = Class, y = DEP_total, fill = Class)) +
  ggplot2::geom_boxplot(alpha = 0.85, outlier.shape = 21, outlier.fill = "white") +
  ggplot2::stat_summary(fun.data = ggplot2::mean_cl_normal, geom = "errorbar",
                        width = 0.2, color = "black", linewidth = 0.8) +
  ggplot2::stat_summary(fun = mean, geom = "point",
                        shape = 23, size = 4, fill = "yellow", color = "black") +
  ggplot2::scale_fill_brewer(palette = "Set1") +
  ggplot2::labs(x = "Latent Class", y = "Depression Total Score (PHQ-9)", fill = "Class") +
  theme_manuscript(base_size = 16) +
  ggplot2::theme(legend.position = "none")

ggplot2::ggsave(
  file.path(FIG_DIR, "Figure3_PHQ9_boxplot_by_class.png"),
  p_dep_box,
  width = 7,
  height = 5,
  dpi = 600
)

## ---- Sensitivity analysis: LPA excluding P ----
lpa_sens_dat <- lpa_dat %>% dplyr::select(id, E, N, L, DEP_total)

set.seed(SEED)
model_list_sens <- tidyLPA::estimate_profiles(
  lpa_sens_dat[, c("E", "N", "L")],
  n_profiles = 1:5,
  models = 1,
  seed = SEED,
  nrep = 5
)

fit_info_sens <- tidyLPA::get_fit(model_list_sens)
save_csv(fit_info_sens, "Sensitivity_LPA_excluding_P_fit_indices.csv")

idx_sens <- which(fit_info_sens$Classes == 4)
model_sens <- model_list_sens[[idx_sens]]
dat_sens <- tidyLPA::get_data(model_sens)
dat_sens$id <- lpa_sens_dat$id

full_dat_sens <- lpa_sens_dat %>%
  dplyr::left_join(dat_sens[, c("id", "Class")], by = "id") %>%
  dplyr::mutate(Class = factor(Class))

desc_sens <- full_dat_sens %>%
  dplyr::group_by(Class) %>%
  dplyr::summarise(
    n = dplyr::n(),
    E_mean = mean(E, na.rm = TRUE),
    N_mean = mean(N, na.rm = TRUE),
    L_mean = mean(L, na.rm = TRUE),
    DEP_mean = mean(DEP_total, na.rm = TRUE),
    .groups = "drop"
  )
save_csv(desc_sens, "Sensitivity_LPA_excluding_P_class_summary.csv")

fit_aov_sens <- stats::aov(DEP_total ~ Class, data = full_dat_sens)
print(summary(fit_aov_sens))
print(stats::TukeyHSD(fit_aov_sens))

saveRDS(
  list(
    dat = dat,
    lpa_dat = lpa_dat,
    model_list = model_list,
    fit_info = fit_info,
    model4 = model4,
    full_dat4 = full_dat4,
    profile_summary = profile_summary,
    table1_panelA = table1_panelA,
    tukey_table = tukey_table,
    full_dat_sens = full_dat_sens,
    fit_info_sens = fit_info_sens
  ),
  file = file.path(RDS_DIR, "lpa_results.rds")
)
