# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(tidyverse)
    library(reshape2)
    library(ComplexHeatmap)
    library(RColorBrewer)
    library(viridis)
    library(ggsignif)
    library(patchwork)
})

PROCDATA_DIR <- "data/procdata/"
RESULTS_DIR <- "data/results/"

source("workflow/scripts/2-AA-analysis/utils/plots.R")
source("workflow/scripts/2-AA-analysis/utils/utils.R")
source("workflow/scripts/utils/palettes.R")
source("workflow/scripts/utils/cohorts.R")

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

# remove unneeded columns from AA_results
AA_results <- AA_results[,-c(6,9,16:28)]

# fix column names
colnames(AA_results) <- gsub(" ", "_", colnames(AA_results))

# -------------------------------------------------------------
# Oncogene-specific analysis
# -------------------------------------------------------------

###########################################################
# Split genes and extract gene features
###########################################################

# get oncogenes
amplified_onco <- AA_results[AA_results$Oncogene != "[]",] # 2139 amplicons with oncogene
oncogenes <- unlist(
  strsplit(gsub("\\[|\\]|'", "", amplified_onco$Oncogenes), ",\\s*")
) |> unique() # 813 unique amplified oncogenes

# create new dataframe
oncogene_df <- data.frame(matrix(nrow=0, ncol=9))
for (oncogene in oncogenes) {
    samples <- amplified_onco[grep(oncogene, amplified_onco$Oncogene),c(1:2, 5, 8:9, 11:12, 14)]
    samples$oncogene <- oncogene
    oncogene_df <- rbind(oncogene_df, samples)
}

###########################################################
# Top amplified oncogenes across all amplicon classes
###########################################################

count_oncogene <- data.frame(table(oncogene_df$oncogene))
count_oncogene <- count_oncogene[order(count_oncogene$Freq, decreasing = TRUE),]

# subset for oncogenes amplified in >= 20 amplicons (n=34)
top_oncogene <- count_oncogene[count_oncogene$Freq >= 20,]

# AOV + tukey test
sig_ecdna_linear <- c()
for (gene in as.character(top_oncogene$Var1)) {
    subset_df <- oncogene_df[oncogene_df$oncogene == gene,]
    # feature median CN AOV and tukey test
    aov.res <- aov(Feature_median_copy_number ~ Classification, data = subset_df)
    tukey <- TukeyHSD(aov.res)$Classification |> as.data.frame()
    sig_tukey <- tukey[tukey$"p adj" < 0.05,]
    if (nrow(sig_tukey) > 0) {
        if ("Linear-ecDNA" %in% rownames(sig_tukey)) {
            sig_ecdna_linear <- c(sig_ecdna_linear, gene)
        }
    }
}

plot_top_amp_oncogenes(top_oncogene, oncogene_df, sig_ecdna_linear)

###########################################################
# Top prevalent oncogenes per amplicon class
###########################################################

plot_top_amp_oncogenes_class(oncogene_df, "BFB")
plot_top_amp_oncogenes_class(oncogene_df, "Complex-non-cyclic")
plot_top_amp_oncogenes_class(oncogene_df, "Linear")
plot_top_amp_oncogenes_class(oncogene_df, "ecDNA")

# -------------------------------------------------------------
# All gene anlaysis
# -------------------------------------------------------------

###########################################################
# Split genes and extract gene features
###########################################################

# get oncogenes
amplified_gene <- AA_results[AA_results$All_genes != "[]",] # 5864 amplicons with oncogene
all_genes <- unlist(
  strsplit(gsub("\\[|\\]|'", "", amplified_gene$All_genes), ",\\s*")
) |> unique() # 16805 unique amplified oncogenes

# create new dataframe
all_genes_df <- data.frame(matrix(nrow=0, ncol=10))
for (gene in all_genes) {
    samples <- amplified_gene[grep(gene, amplified_gene$All_genes),c(1:2, 5, 8:9, 11:12, 14)]
    samples$gene <- gene
    samples$oncogene <- ifelse(gene %in% oncogenes, "Oncogene", "Not")
    all_genes_df <- rbind(all_genes_df, samples)
}

saveRDS(all_genes_df, file = paste0(PROCDATA_DIR, "AA_exploration/all_genes_df.rds"))

###########################################################
# Top amplified genes across all amplicon classes
###########################################################

count_gene <- data.frame(table(all_genes_df$gene))
count_gene <- count_gene[order(count_gene$Freq, decreasing = TRUE),]
count_gene$oncogene <- ifelse(count_gene$Var1 %in% oncogenes, "Oncogene", "Not")

# subset for top 50
top_gene <- count_gene[1:50,]

# AOV + tukey test
sig_ecdna_linear <- c()
for (gene in as.character(top_gene$Var1)) {
    subset_df <- all_genes_df[all_genes_df$gene == gene,]
    # feature median CN AOV and tukey test
    aov.res <- aov(Feature_median_copy_number ~ Classification, data = subset_df)
    tukey <- TukeyHSD(aov.res)$Classification |> as.data.frame()
    sig_tukey <- tukey[tukey$"p adj" < 0.05,]
    if (nrow(sig_tukey) > 0) {
        if ("Linear-ecDNA" %in% rownames(sig_tukey)) {
            sig_ecdna_linear <- c(sig_ecdna_linear, gene)
        }
    }
} # 0 sig results here which is weird bc MYC is here for example

plot_top_amp_genes(top_gene, all_genes_df)

###########################################################
# Top prevalent oncogenes per amplicon class
###########################################################

plot_top_amp_all_genes_class(all_genes_df, "BFB")
plot_top_amp_all_genes_class(all_genes_df, "Complex-non-cyclic")
plot_top_amp_all_genes_class(all_genes_df, "Linear")
plot_top_amp_all_genes_class(all_genes_df, "ecDNA")