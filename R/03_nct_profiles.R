## =========================================================
## 03_nct_profiles.R
## Network Comparison Tests across latent profiles
## =========================================================

source(file.path("R", "00_setup.R"))

lpa_results <- readRDS(file.path(RDS_DIR, "lpa_results.rds"))
dat <- lpa_results$dat
full_dat4 <- lpa_results$full_dat4

dat_phq <- dat %>% dplyr::select(id, dplyr::all_of(PHQ_VARS))
full_dat_phq <- full_dat4 %>% dplyr::left_join(dat_phq, by = "id")

dat_c1 <- full_dat_phq %>% dplyr::filter(as.character(Class) == "1") %>% dplyr::select(dplyr::all_of(PHQ_VARS))
dat_c2 <- full_dat_phq %>% dplyr::filter(as.character(Class) == "2") %>% dplyr::select(dplyr::all_of(PHQ_VARS))
dat_c3 <- full_dat_phq %>% dplyr::filter(as.character(Class) == "3") %>% dplyr::select(dplyr::all_of(PHQ_VARS))
dat_c4 <- full_dat_phq %>% dplyr::filter(as.character(Class) == "4") %>% dplyr::select(dplyr::all_of(PHQ_VARS))

set.seed(SEED)
NCT_ITER <- 1000

nct_3v1 <- NetworkComparisonTest::NCT(dat_c3, dat_c1, it = NCT_ITER, test.edges = FALSE, test.centrality = FALSE)
nct_3v2 <- NetworkComparisonTest::NCT(dat_c3, dat_c2, it = NCT_ITER, test.edges = FALSE, test.centrality = FALSE)
nct_3v4 <- NetworkComparisonTest::NCT(dat_c3, dat_c4, it = NCT_ITER, test.edges = FALSE, test.centrality = FALSE)

extract_nct <- function(nct_obj, comparison) {
  data.frame(
    Comparison = comparison,
    M_statistic = nct_obj$nwinv.real,
    M_p = nct_obj$nwinv.pval,
    Global_strength_group1 = nct_obj$glstrinv.real[1],
    Global_strength_group2 = nct_obj$glstrinv.real[2],
    S_statistic = abs(diff(nct_obj$glstrinv.real)),
    S_p = nct_obj$glstrinv.pval,
    stringsAsFactors = FALSE
  )
}

nct_summary <- dplyr::bind_rows(
  extract_nct(nct_3v1, "Class 3 vs Class 1"),
  extract_nct(nct_3v2, "Class 3 vs Class 2"),
  extract_nct(nct_3v4, "Class 3 vs Class 4")
)

print(nct_summary)
save_csv(nct_summary, "NCT_summary_Class3_vs_other_profiles.csv")

saveRDS(
  list(nct_3v1 = nct_3v1, nct_3v2 = nct_3v2, nct_3v4 = nct_3v4, nct_summary = nct_summary),
  file = file.path(RDS_DIR, "nct_results.rds")
)
