# Run first three steps from the repository root.
# I keep this small on purpose: it is only for checking the early data flow.

source("scripts/01_prepare_TCGA_data.R")
source("scripts/02_ESTIMATE_score.R")
source("scripts/03_survival_and_clinical.R")
