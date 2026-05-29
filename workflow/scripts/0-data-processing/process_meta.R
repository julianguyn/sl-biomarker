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
    Patient_ID = sub("^([A-Za-z0-9]+_[0-9]+)_.*$", "\\1", all_samples$sample),
    PatCan_ID = sub("^([A-Za-z0-9]+_[0-9]+_[A-Za-z]+_[A-Z])_.*$", "\\1", all_samples$sample),
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
# Create metadata of available RSEM
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
rsem_meta$PatCan_ID <- sub("^([A-Za-z0-9]+_[0-9]+_[A-Za-z]+_[A-Z])_.*$", "\\1", rsem_meta$Omics_ID)

# 1690 unique sample IDs from PM_meta
# 461 from batch1, 1118 from batch2 = 1579 total

table(rsem_meta$Patient_ID %in% PM_meta$Patient_ID) #all 1579 match

###########################################################
# Check duplicates
###########################################################

# one duplicate in the PM_meta patient IDs (P and M)
duplicated_PM_id <- PM_meta$Patient_ID[duplicated(PM_meta$Patient_ID)]

# 58 duplicates in the RSEM
duplicated_PM_rna <- rsem_meta$Patient_ID[duplicated(rsem_meta$Patient_ID)]

# ------ check duplicates
tt <- rsem_meta[rsem_meta$Patient_ID %in% duplicated_PM_rna,]
tt <- tt[order(tt$Patient_ID),]

b1 <- tt[tt$Batch == "Batch1",] #45
b2 <- tt[tt$Batch == "Batch2",] #70

# 45 are duplicates in the batches (they have the same Omics_ID)
two_batches <- duplicated_PM_rna[which(duplicated_PM_rna %in% b1$Patient_ID & duplicated_PM_rna %in% b2$Patient_ID)]
table(b1$Omics_ID[b1$Patient_ID %in% two_batches] == b2$Omics_ID[b2$Patient_ID %in% two_batches])

# remove these duplicates
rsem_meta <- rsem_meta[-which(rsem_meta$Patient_ID %in% two_batches & rsem_meta$Batch == "Batch2"),]
# 1534 left in rsem_meta

# remaining duplicates are all in b2
b2 <- b2[-which(b2$Patient_ID %in% two_batches),] #25, 12 unique PatientIDs

###########################################################
# Try mapping to the ecDNA IDs
###########################################################

rsem_meta$match_WG <- sub("_WT", "_WG", sub("RNA", "DNA", rsem_meta$Omics_ID)) # these are all unique
table(rsem_meta$match_WG %in% PM_meta$Sample_ID) # 1472 match, 62 don't match

###### check the ones that don't match
no_match <- rsem_meta[-which(rsem_meta$match_WG %in% PM_meta$Sample_ID),] #62
dup <- no_match$Patient_ID[duplicated(no_match$Patient_ID)] # 1 duplicate

table(no_match$Patient_ID %in% PM_meta$Patient_ID) # all true lol

no_match$Matched_ID <- NA
for (sample in no_match$Patient_ID) {

    match <- PM_meta[PM_meta$Patient_ID == sample,]
    if (nrow(match) == 1) {
        no_match$Matched_ID[no_match$Patient_ID == sample] <- match$Sample_ID
    } else {
        cat("----- Sample:", sample, "\n", nrow(match))
    }
}

#### check the ones that do match
matched <- rsem_meta[rsem_meta$match_WG %in% PM_meta$Sample_ID,] #1472

# check if any of these had an indirect match wtih a B2 dup
to_check <- no_match[no_match$Patient_ID %in% unique(b2$Patient_ID),] # 13 cases here
to_check[to_check$Matched_ID %in% matched$match_WG,] |> dim() # all 13 had a better direct match

b2$removed <- ifelse(b2$Omics_ID %in% to_check$Omics_ID, "Removed", "Kept")
no_match <- no_match[-which(no_match$Omics_ID %in% to_check$Omics_ID),] # 49 manually mapped remaining

tt <- no_match[,c(2,5,6)] 
write.csv(tt, file = "metadata/manual_rna_ecdna_match.csv", quote = FALSE, row.names = FALSE)

# merge back (1521 samples matched to RNA)
length(unique(matched$Omics_ID))
length(matched$Omics_ID)

length(unique(no_match$Omics_ID))
length(no_match$Omics_ID)

length(unique(c(matched$Omics_ID, no_match$Omics_ID)))

###########################################################
# Clean up RSEM Meta
###########################################################s
 
# 1521 matches (1472 direct, 49 indirect), removed 13 duplicates
rsem_meta$Status <- ifelse(
    rsem_meta$Omics_ID %in% matched$Omics_ID, "Direct Match",
    ifelse(
        rsem_meta$Omics_ID %in% to_check$Omics_ID, "Removed (Dup)", "Indirect Match"
    )
)

for (i in 1:nrow(rsem_meta)) {
    sample <- rsem_meta$Omics_ID[i]
    if (rsem_meta$Status[i] == "Indirect Match") {
        rsem_meta$match_WG[i] <- no_match$Matched_ID[no_match$Omics_ID == sample]
    } else if (rsem_meta$Status[i] == "Removed (Dup)") {
        rsem_meta$match_WG[i] <- "Removed"
    }
}
table(rsem_meta$match_WG %in% PM_meta$Sample_ID) # 1521 match, 13 dups

###########################################################
# Add RNA info to PM_meta
###########################################################s
 
PM_meta$RNA <- "missing"
for (sample in PM_meta$Sample_ID) {
    if (sample %in% rsem_meta$match_WG) {
        PM_meta$RNA[PM_meta$Sample_ID == sample] <- rsem_meta$Omics_ID[rsem_meta$match_WG == sample]
    }
}
dim(PM_meta[PM_meta$RNA == "missing",]) #169 cases with missing RNA

write.csv(PM_meta, file = "metadata/2026-05-29_PM_meta.csv", quote = FALSE, row.names = FALSE)

missing_RNA <- PM_meta[PM_meta$RNA == "missing",]
write.csv(missing_RNA, file = "metadata/2026-05-29_PM_meta_missing_RNA.csv", quote = FALSE, row.names = FALSE)

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