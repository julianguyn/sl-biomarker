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
omics <- read.table("metadata/meta.tsv", header = TRUE)

###########################################################
# Create full metadata
###########################################################

meta <- data.frame(
    Centre = all_samples$Centre,
    Cohort = all_samples$cohort,
    Sample_ID = all_samples$sample,
    MOHCCN_ID = NA,
    Patient_ID = sub("^([A-Za-z0-9]+_[0-9]+)_.*$", "\\1", all_samples$sample)
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
# Make life easier by splitting PM and BC
###########################################################

PM_meta <- meta[meta$Centre == "PM2C",]
BC_meta <- meta[meta$Centre == "BCCA",]

###########################################################
# Map RNA to PM
###########################################################

load("metadata/PM2C_RSEM.RData")
# included here is batch1 and batch2, two vectors of strings

rsem_meta <- data.frame(
    Batch = c(
        rep("Batch1", length(batch1)),
        rep("Batch2", length(batch2))
    ),
    Omics_ID = sub("\\.RSEM.*", "", c(batch1, batch2))
)
rsem_meta$Patient_ID <- sub("^([A-Za-z0-9]+_[0-9]+)_.*$", "\\1", rsem_meta$Omics_ID)


# 1690 unique sample IDs from PM_meta
# 451 from batch1, 1118 from batch2 = 1569 total

# trying to map the IDs 

table(rsem_meta$Omics_ID %in% PM_meta$Sample_ID)


table(sub("_WT_", "_WG_", batch1) %in% PM_meta$Sample_ID)

###########################################################
# Mapping omics to ecDNA
###########################################################

# first pass mapping
meta$RNA <- "missing"
meta$MUT <- "missing"
meta$Omics_ID <- NA
for (sample in intersect(meta$Sample_ID, omics$OICR_ID)) {
    meta$MOHCCN_ID[meta$Sample_ID == sample] <- omics$PM2C_SAMPLE_ID[omics$OICR_ID == sample]
    meta$RNA[meta$Sample_ID == sample] <- omics$RNA[omics$OICR_ID == sample]
    meta$MUT[meta$Sample_ID == sample] <- omics$MUT[omics$OICR_ID == sample]
    meta$Omics_ID[meta$Sample_ID == sample] <- sample
}

# helper function for manual pass
map_IDs <- function(meta, id_in_meta, id_in_omics) {
    meta$MOHCCN_ID[meta$Sample_ID == id_in_meta] <- omics$PM2C_SAMPLE_ID[omics$OICR_ID == id_in_omics]
    meta$RNA[meta$Sample_ID == id_in_meta] <- omics$RNA[omics$OICR_ID == id_in_omics]
    meta$MUT[meta$Sample_ID == id_in_meta] <- omics$MUT[omics$OICR_ID == id_in_omics]
    meta$Omics_ID[meta$Sample_ID == id_in_meta] <- id_in_omics
    return(meta)
}

# second pass with manual mapping
for (sample in names(ID_MAPPING)) {
    meta <- map_IDs(meta, sample, ID_MAPPING[[sample]])
}

###########################################################
# Load in BC Mapping
###########################################################

bc_meta <- fread("metadata/MOHCCN-BC_sample_file_map_merged.tsv", data.table = FALSE)
load("metadata/RNASeq.RData")

# get IDs shared between the RNA data and the BC metadata
bc_meta$have_RNA <- ifelse(bc_meta$library %in% ids, "RNA", "no RNA")
have_RNA <- bc_meta[bc_meta$have_RNA == "RNA",]

# get IDs shared between the ecDNA data and RNA data
bc_ids <- intersect(have_RNA$donor_study_id, meta$Sample_ID)

for (sample in bc_ids) {
    library_id <- have_RNA$library[have_RNA$donor_study_id == sample]
    if (length(library_id) > 1) library_id <- paste0(library_id[1], "_", library_id[2])
    meta$Omics_ID[meta$Sample_ID == sample] <- library_id
}

###########################################################
# Save metadata
###########################################################

write.table(meta, file = "metadata/2026-05-28_meta_ecDNA.tsv", quote = FALSE, row.names = FALSE)

bc_no_ecDNA <- have_RNA[-which(have_RNA$donor_study_id %in% meta$Sample_ID),]
save(bc_no_ecDNA, removed_duplicates, file = "metadata/sample_ids_removed.RData")