# Personality Profiles and Depressive Symptom Networks

This repository contains the R analysis scripts used for the manuscript:

**From Personality Profiles to Depressive Symptom Networks: Exploratory In Silico Projections of Potential Symptom Targets among Emerging Adults**

The repository is intended to support analytic transparency and reproducibility. 
## Repository structure

```text
R/
  00_setup.R
  01_lpa_profiles.R
  02_network_highrisk.R
  03_nct_profiles.R
  04_insilico_simulation.R
data/
  .gitkeep
outputs/
  tables/
  figures/
  rds/
```

## Running the scripts

Run the scripts from the repository root in the following order:

```r
source("R/01_lpa_profiles.R")
source("R/02_network_highrisk.R")
source("R/03_nct_profiles.R")
source("R/04_insilico_simulation.R")
```

The file `R/00_setup.R` is sourced automatically by the analysis scripts and contains shared settings, package loading, variable labels, and helper functions.

## Data

To run the scripts, place the cleaned analysis dataset at:

```text
data/analysis_data.sav
```

or set the data path manually before running the scripts:

```r
Sys.setenv(STUDY_DATA_PATH = "path/to/your/analysis_data.sav")
```

The scripts assume that missing data have already been handled in the cleaned analysis dataset, as described in the manuscript.

## Required variables

The scripts are designed for a de-identified, cleaned analysis dataset. The dataset is not included in this repository. To run the scripts, users should place the analysis dataset in the data/ folder and ensure that it contains the following variables:

E: Extraversion score
N: Neuroticism score
P: Psychoticism score
L: Lie/Social Desirability score
DEP_total: PHQ-9 total score
D1–D9: PHQ-9 item scores
## Main analysis settings

- Random seed: `2025`
- LPA: `tidyLPA` Model 1; one- to five-profile solutions
- Retained profile solution: four profiles
- Main network: EBICglasso with tuning parameter `gamma = 0.50`
- Network Comparison Tests: `1,000` permutations
- Bootstrap analyses: `2,000` bootstrap samples
- Simulation perturbation: hypothetical 1-SD reduction, truncated to the PHQ-9 item range of 0-3

## Outputs

Generated files are saved locally to:

```text
outputs/tables/
outputs/figures/
outputs/rds/
```

The `.gitignore` file excludes raw data and generated outputs from version control. Do not commit the dataset or locally generated result files.
