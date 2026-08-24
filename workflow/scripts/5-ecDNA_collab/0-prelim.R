##### quickly check associations

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

# load in RNA
PM_rna <- readRDS("data/procdata/RNA_TPM/PM_RNA.rds")

# load in AA results
AA_results <- readRDS(paste0(PROCDATA_DIR, "AA_output/AA_results.rds"))

# load in all samples
all_samples <- readRDS(paste0(PROCDATA_DIR, "AA_output/sample_df.rds"))

# process cohort names (from utils/get_tables.R)
AA_results <- process_cohort_names(AA_results, res = TRUE)
all_samples <- process_cohort_names(all_samples)

###########################################################
# Preprocessing
###########################################################

# remove duplicates
all_samples <- all_samples[-which(all_samples$duplicated == "duplicated_removed"),]

# remove missing amplicons
all_samples <- all_samples[-which(all_samples$result_table == "missing"),]

# remove NA amplicons
AA_results <- AA_results[!is.na(AA_results$'AA amplicon number'),]
all_samples <- all_samples[all_samples$sample %in% AA_results$'Sample name',]

# add patient ID column
AA_results$Patient_ID <- sub("^([A-Za-z0-9]+_[0-9]+)_.*$", "\\1", AA_results$'Sample name')

# get ecDNA
ecDNA <- AA_results[
    AA_results$Classification == "ecDNA",
    c("Centre", "cohort", "Sample name", "Patient_ID", "Feature ID", "Oncogenes", "Feature median copy number", "Feature maximum copy number")
] |> as.data.frame()

###########################################################
# Quick numbers
###########################################################

# number of patients with ecDNA
dim(ecDNA)
length(unique(ecDNA$Patient_ID))
myc <- ecDNA[grep("\'MYC\'", ecDNA$Oncogenes),]
dim(myc)
length(unique(myc$Patient_ID))


# figure out cancer type distribution
type <- sub("_.*", "", sub("^([A-Za-z0-9]+_[0-9]+)_", "", PM_meta$PatCan_ID))

ecDNA_type <- sub("^([A-Za-z0-9]+_[0-9]+_[A-Za-z]+)_.*$", "\\1", ecDNA$'Sample name')
ecDNA_type <- sub("_.*", "", sub("^([A-Za-z0-9]+_[0-9]+)_", "", ecDNA_type))

myc_meta <- PM_meta[PM_meta$Patient_ID %in% myc$Patient_ID,]

# --- remove later, keep just PM cases for now
ecDNA <- ecDNA[ecDNA$Centre == "PM2C",] #700

###########################################################
# Correlation: PVT1 expression and ecDNA copy number
###########################################################

# keep PVT1 expression
pvt1 <- "ENSG00000249859.11"
pvt1_rna <- as.data.frame(t(PM_rna[pvt1,]))
pvt1_rna$Patient_ID <- sub("^([A-Za-z0-9]+_[0-9]+)_.*$", "\\1", rownames(pvt1_rna))


### ---- with PVT1 ecDNA

# ecDNA with PVT1
pvt1_ecDNA <- ecDNA$Patient_ID[grep("PVT1", ecDNA$Oncogenes)] #25 patients


pvt1_rna$have_ecDNA <- ifelse(
    pvt1_rna$Patient_ID %in% ecDNA$Patient_ID, 
    ifelse(pvt1_rna$Patient_ID %in% pvt1_ecDNA, "PVT1_ecDNA+", "PVT1_ecDNA-"), 
    "ecDNA-"
)

# PVT1 expression in patients with and without ecDNA (ecDNA with or without PVT1)
ggplot(pvt1_rna, aes(x = have_ecDNA, y = log2(ENSG00000249859.11 + 1))) +
    geom_boxplot() + geom_jitter(width = 0.2)

# PVT1 ecDNA - expr correlation with ecDNA cn
pvt1_ecDNA_cn <- ecDNA[grep("PVT1", ecDNA$Oncogenes),]
pvt1_ecDNA_cn$RNA <- pvt1_rna$ENSG00000249859.11[match(pvt1_ecDNA_cn$Patient_ID, pvt1_rna$Patient_ID)]


ggplot(pvt1_ecDNA_cn, aes(x = RNA, y = .data[['Feature median copy number']])) +
    geom_point()

ggplot(pvt1_ecDNA_cn, aes(x = RNA, y = .data[['Feature maximum copy number']])) +
    geom_point()


### ---- with MYC ecDNA

# ecDNA with PVT1
myc_ecDNA <- ecDNA$Patient_ID[grep("\'MYC\'", ecDNA$Oncogenes)] #25 patients


pvt1_rna$have_MYC_ecDNA <- ifelse(
    pvt1_rna$Patient_ID %in% ecDNA$Patient_ID, 
    ifelse(pvt1_rna$Patient_ID %in% myc_ecDNA, "MYC_ecDNA+", "MYC_ecDNA-"), 
    "ecDNA-"
)

# PVT1 expression in patients with and without ecDNA (ecDNA with or without PVT1)
ggplot(pvt1_rna, aes(x = have_MYC_ecDNA, y = log2(ENSG00000249859.11 + 1))) +
    geom_boxplot() + geom_jitter(width = 0.2)

# PVT1 ecDNA - expr correlation with ecDNA cn
myc_ecDNA_cn <- ecDNA[grep("\'MYC\'", ecDNA$Oncogenes),]
myc_ecDNA_cn$RNA <- pvt1_rna$ENSG00000249859.11[match(myc_ecDNA_cn$Patient_ID, pvt1_rna$Patient_ID)]


ggplot(myc_ecDNA_cn, aes(x = RNA, y = .data[['Feature median copy number']])) +
    geom_point()

ggplot(pvt1_ecDNA_cn, aes(x = RNA, y = .data[['Feature maximum copy number']])) +
    geom_point()
