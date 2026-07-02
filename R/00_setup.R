## =========================================================
## 00_setup.R
## Shared setup
## =========================================================

required_packages <- c(
  "haven", "dplyr", "tidyr", "ggplot2", "tidyLPA",
  "qgraph", "bootnet", "NetworkComparisonTest",
  "effectsize", "igraph"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Please install required packages: ", paste(missing_packages, collapse = ", "))
}

library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyLPA)
library(qgraph)
library(bootnet)
library(NetworkComparisonTest)
library(effectsize)
library(igraph)

SEED <- 2025
set.seed(SEED)

DATA_PATH <- Sys.getenv("STUDY_DATA_PATH", unset = file.path("data", "analysis_data.sav"))

TABLE_DIR <- file.path("outputs", "tables")
FIG_DIR <- file.path("outputs", "figures")
RDS_DIR <- file.path("outputs", "rds")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

EPQ_COLS <- 130:133
DEP_COL <- 134
PHQ_VARS <- paste0("D", 1:9)
HIGH_RISK_CLASS <- "3"

CLASS_LABELS <- data.frame(
  Class = as.character(1:4),
  Profile_label = c(
    "Mildly Dysregulated",
    "Adaptive-Conforming",
    "Vulnerable-Detached",
    "Stable-Detached"
  ),
  stringsAsFactors = FALSE
)

PHQ_LABELS <- c(
  D1 = "Anhedonia",
  D2 = "Depressed mood",
  D3 = "Sleep problems",
  D4 = "Fatigue",
  D5 = "Appetite change",
  D6 = "Low self-esteem",
  D7 = "Concentration problems",
  D8 = "Psychomotor change",
  D9 = "Suicidal ideation"
)

PHQ_AXIS_LABELS <- c(
  D1 = "D1 Anhedonia",
  D2 = "D2 Depressed mood",
  D3 = "D3 Sleep problems",
  D4 = "D4 Fatigue",
  D5 = "D5 Appetite change",
  D6 = "D6 Low self-esteem",
  D7 = "D7 Concentration",
  D8 = "D8 Psychomotor change",
  D9 = "D9 Suicidal ideation"
)

load_clean_data <- function(data_path = DATA_PATH) {
  if (!file.exists(data_path)) {
    stop("Data file not found: ", data_path,
         "\nPlace the cleaned .sav file at data/analysis_data.sav or set STUDY_DATA_PATH.")
  }

  dat <- haven::read_sav(data_path)
  names(dat)[EPQ_COLS] <- c("E", "N", "P", "L")
  names(dat)[DEP_COL] <- "DEP_total"

  dat <- dat %>% dplyr::mutate(id = dplyr::row_number())

  required_vars <- c("id", "E", "N", "P", "L", "DEP_total", PHQ_VARS)
  missing_vars <- setdiff(required_vars, names(dat))
  if (length(missing_vars) > 0) {
    stop("Missing required variables: ", paste(missing_vars, collapse = ", "))
  }

  dat
}

make_lpa_data <- function(dat) {
  lpa_dat <- dat %>% dplyr::select(id, E, N, P, L, DEP_total)
  if (anyNA(lpa_dat)) {
    stop("Missing values detected in LPA variables. Use the cleaned dataset after missing-data handling.")
  }
  lpa_dat
}

save_csv <- function(x, file) {
  utils::write.csv(x, file.path(TABLE_DIR, file), row.names = FALSE)
}

theme_manuscript <- function(base_size = 14) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.title = ggplot2::element_text(face = "bold")
    )
}
