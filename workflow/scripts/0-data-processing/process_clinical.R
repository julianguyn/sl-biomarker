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

# load in primary diagnosis information
primary <- read.csv("data/rawdata/clinical_data_20260923/primary_diagnoses.csv")

# load in ID mapping file
map <- read.table("data/rawdata/ID_mapping_20260925/MOHCCN_PM2C_ID_map.tsv", header = T)
map$PM2C_match_ID <- sub("^(([^-]+-){2}[^-]+).*$", "\\1", map$PM2C_SAMPLE_ID)

# map IDs
primary$sample <- primary$submitter_donor_id
for (i in 1:nrow(primary)) {
    sample <- primary$submitter_donor_id[i]
    if (sample %in% map$PM2C_match_ID) {
        id <- map$OICR_SAMPLE_FILE[map$PM2C_match_ID == sample]
        primary$sample[i] <- id
    }
}

###########################################################
# Overlap sampleIDs in primary
###########################################################

primary <- primary[,c("program_id", "submitter_donor_id", "sample", "primary_site")]
search <- all_samples[,c("cohort", "sample", "Centre")]
search <- search[-which(search$sample == "PM2C_Batch_2"),]
search$cohort <- as.character(search$cohort)

search$primary <- ifelse(search$sample %in% primary$sample, primary$primary_site, "missing")

###########################################################
# Get missing data
###########################################################

missing <- search[search$primary == "missing",]

# get the missing data
PM_missing <- missing[missing$Centre == "PM2C",]
BC_missing <- missing[missing$Centre == "BCCA",]

###########################################################
# Second round of matching
###########################################################

PM_missing$patientID <- sub("^((?:[^_]+_){1}[^_]+).*", "\\1", PM_missing$sample)
map <- map[map$OICR_sID %in% PM_missing$patientID,]

tt <- primary[primary$sample %in% map$OICR_SAMPLE_FILE,]

to_map <- map[map$OICR_SAMPLE_FILE %in% primary$sample,]
to_map$primary <- primary$primary_site[match(to_map$OICR_SAMPLE_FILE, primary$sample)]

# map primary
PM_found <- PM_missing[PM_missing$patientID %in% to_map$OICR_sID,]
PM_found$primary <- to_map$primary[match(PM_found$patientID, to_map$OICR_sID)]

# get samples still missing
PM_missing <- PM_missing[-which(PM_missing$patientID %in% PM_found$patientID),]

###########################################################
# Save
###########################################################

#write.csv(search, file = "data/procdata/clinical/coded_missing_primary.csv", quote = FALSE, row.names = FALSE)
#write.csv(missing, file = "data/procdata/clinical/missing_primary.csv", quote = FALSE, row.names = FALSE)
