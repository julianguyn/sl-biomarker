# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

PROCDATA_DIR <- "data/procdata/"
RESULTS_DIR <- "data/results/"

source("workflow/scripts/utils/palettes.R")
source("workflow/scripts/utils/cohorts.R")
source("workflow/scripts/2-AA-analysis/utils/utils.R")

###########################################################
# Load in data
###########################################################

# load in metadata
PM_meta <- read.csv("metadata/2026-05-29_PM_meta.csv")
BC_meta <- read.csv("metadata/2026-05-29_BC_meta.csv")

# load in AA results
AA_results <- readRDS(paste0(PROCDATA_DIR, "AA_output/AA_results.rds"))

# load in all samples
all_samples <- readRDS(paste0(PROCDATA_DIR, "AA_output/sample_df.rds"))

# process cohort names (from utils/get_tables.R)
AA_results <- process_cohort_names(AA_results, res = TRUE)
all_samples <- process_cohort_names(all_samples)

###########################################################
# Load in clinical data
###########################################################

primary <- read.csv("data/rawdata/clinical_data_20260923/primary_diagnoses.csv")

###########################################################
# Overlap sampleIDs in primary
###########################################################

primary <- primary[,c("program_id", "submitter_donor_id", "primary_site")]
search <- all_samples[,c("cohort", "sample", "Centre")]
search$cohort <- as.character(search$cohort)

search$primary <- ifelse(search$sample %in% primary$submitter_donor_id, search$sample, "missing")

have <- search[search$primary != "missing",]
missing <- search[search$primary == "missing",]
table(primary$submitter_donor_id %in% search$sample)

write.csv(search, file = "data/procdata/clinical/coded_missing_primary.csv", quote = FALSE, row.names = FALSE)
write.csv(missing, file = "data/procdata/clinical/missing_primary.csv", quote = FALSE, row.names = FALSE)


search_AA <- AA_results[,c("Sample name", "Centre"),]
search_AA$primary <- ifelse(search_AA$'Sample name' %in% primary$submitter_donor_id, search$sample, "missing")
have <- search_AA[search_AA$primary != "missing",]
missing <- search_AA[search_AA$primary == "missing",]
