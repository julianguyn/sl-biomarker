##### quickly check associations

# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(dplyr)
    library(purrr)
    library(broom)
    library(ggrepel)
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

# load in mutations
PM_mut <- read.table("data/procdata/mutations/PM_mutations.tsv")
rownames(PM_mut) <- gsub("\\.", "-", rownames(PM_mut))

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
# Set up correlation analysis
###########################################################

pvt1_gene <- "ENSG00000249859.11"
tp53_gene <- "ENSG00000141510.17"

# function to get correlation
get_corr <- function(df, label) {
    corr <- cor.test(df$TP53, df$PVT1, method = "spearman")
    to_bind <- data.frame(
        cor = corr$estimate,
        pval = corr$p.value,
        n = nrow(df),
        label = label
    )
    return(to_bind)
}

PM_rna <- as.data.frame(t(PM_rna[rownames(PM_rna) %in% c(pvt1_gene, tp53_gene),]))
colnames(PM_rna) <- c("TP53", "PVT1")

# make dataframe to score results
corr_res <- data.frame(matrix(nrow=0, ncol=4))


###########################################################
# Correlation of all samples
###########################################################

corr_res <- rbind(corr_res, get_corr(PM_rna, "All Samples"))

###########################################################
# Correlation stratified by TP53 mutation status
###########################################################

# get TP53 mutation status
PM_rna$match_id <- sub("_WG", "", rownames(PM_rna))
PM_rna <- PM_rna[PM_rna$match_id %in% rownames(PM_mut),]
PM_rna$TP53_status <- PM_mut$TP53[match(PM_rna$match_id, rownames(PM_mut))]
PM_rna$TP53_status <- ifelse(PM_rna$TP53_status > 0, "Mut", "Wt")

mut <- PM_rna[PM_rna$TP53_status == "Mut",]
wt <- PM_rna[PM_rna$TP53_status == "Wt",]

# get correlations by TP53 mutation status
corr_res <- rbind(corr_res, get_corr(mut, "TP53 Mut"))
corr_res <- rbind(corr_res, get_corr(wt, "TP53 WT"))

###########################################################
# Correlation stratified by TP53 mutation status and MYC ecDNA
###########################################################

# get MYC ecDNA
myc_ecDNA <- ecDNA[grep("MYC", ecDNA$Oncogenes),]
myc_ecDNA <- myc_ecDNA$'Sample name'[-grep("MYCN", myc_ecDNA$Oncogenes)]
myc_ecDNA <- sub("_WG", "", myc_ecDNA)
PM_rna$MYC_ecDNA <- ifelse(PM_rna$match_id %in% myc_ecDNA, "MYC ecDNA", "non-MYC ecDNA")

mut_myc <- PM_rna[PM_rna$TP53_status == "Mut" & PM_rna$MYC_ecDNA == "MYC ecDNA",]
mut_non_myc <- PM_rna[PM_rna$TP53_status == "Mut" & PM_rna$MYC_ecDNA == "non-MYC ecDNA",]
wt_myc <- PM_rna[PM_rna$TP53_status == "Wt" & PM_rna$MYC_ecDNA == "MYC ecDNA",]
wt_non_myc <- PM_rna[PM_rna$TP53_status == "Wt" & PM_rna$MYC_ecDNA == "non-MYC ecDNA",]

# get correlations by TP53 mutation status
corr_res <- rbind(corr_res, get_corr(mut_myc, "TP53 Mut (MYC ecDNA)"))
corr_res <- rbind(corr_res, get_corr(wt_myc, "TP53 WT (MYC ecDNA)"))
corr_res <- rbind(corr_res, get_corr(mut_non_myc, "TP53 Mut (non-MYC ecDNA)"))
corr_res <- rbind(corr_res, get_corr(wt_non_myc, "TP53 WT (non-MYC ecDNA)"))

###########################################################
# Get final dataframe
###########################################################

df <- corr_res[,c("label", "cor", "pval", "n")]
rownames(df) <- NULL
write.csv(df, file = "data/results/data/ecDNA_collab/pvt1_tp53_corr.csv", quote = FALSE, row.names = FALSE)
