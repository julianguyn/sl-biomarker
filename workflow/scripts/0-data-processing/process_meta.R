# compile all available sample data

library(data.table)


PROCDATA_DIR <- "data/procdata/"

source("workflow/scripts/2-AA-analysis/utils/utils.R")
source("workflow/scripts/utils/cohorts.R")
source("workflow/scripts/utils/ID_MAPPING.R")

###########################################################
# Load in data
###########################################################

# load in AA results
AA_results <- readRDS(paste0(PROCDATA_DIR, "AA_output/AA_results.rds"))

# load in all samples
all_samples <- readRDS(paste0(PROCDATA_DIR, "AA_output/sample_df.rds"))

# process cohort names (from utils/get_tables.R)
AA_results <- process_cohort_names(AA_results, res = TRUE)
all_samples <- process_cohort_names(all_samples)

# omics metadata

###########################################################
# Create full metadata
###########################################################

meta <- data.frame(
    Centre = all_samples$Centre,
    Cohort = all_samples$cohort,
    Sample_ID = all_samples$sample,
    MOHCCN_ID = NA,
    Batch = all_samples$batch,
    Duplicated = all_samples$duplicated,
    Results_Table = all_samples$result_table
)
meta$Results_Table <- ifelse(meta$Results_Table == "missing", "Missing", "Available")


# NA amplicons
na_amps <- AA_results$'Sample name'[is.na(AA_results$'AA amplicon number')]
meta$NA_Amp <- ifelse(meta$Sample_ID %in% na_amps, "NA Amplicon Number", "Has Amplicon")

# split duplicates
removed_duplicates <- meta[meta$Duplicated == "duplicated_removed",]
meta <- meta[-which(meta$Duplicated == "duplicated_removed"),]

###########################################################
# Mapping omics to ecDNA
###########################################################

# first pass mapping
meta$RNA <- "missing"
meta$MUT <- "missing"
meta$Omics_ID <- NA
for (sample in intersect(meta$Sample_ID, omics$OICR_ID)) {
    meta$RNA[meta$Sample_ID == sample] <- omics$RNA[omics$OICR_ID == sample]
    meta$MUT[meta$Sample_ID == sample] <- omics$MUT[omics$OICR_ID == sample]
    meta$Omics_ID[meta$Sample_ID == sample] <- sample
}

# helper function for manual pass
map_IDs <- function(meta, id_in_meta, id_in_omics) {
    meta$RNA[meta$Sample_ID == id_in_meta] <- omics$RNA[omics$OICR_ID == id_in_omics]
    meta$MUT[meta$Sample_ID == id_in_meta] <- omics$MUT[omics$OICR_ID == id_in_omics]
    meta$Omics_ID[meta$Sample_ID == id_in_meta] <- id_in_omics
    return(meta)
}

# second pass with manual mapping
for (sample in names(ID_MAPPING)) {
    meta <- map_IDs(meta, sample, ID_MAPPING[[sample]])
}
