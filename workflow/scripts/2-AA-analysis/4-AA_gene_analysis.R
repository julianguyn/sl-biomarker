# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggridges)
    library(dplyr)
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

# ----------------------------------------------------------------------
# Gene Exploration
# ----------------------------------------------------------------------

###########################################################
# Number of genes in each sample
###########################################################

# get number of amplified oncogenes
AA_results$n_genes <- lengths(
    strsplit(gsub("\\[|\\]|'", "", AA_results$Oncogenes), ",\\s*")
)

# get number of amplified genes
AA_results$n_genes_all <- lengths(
    strsplit(gsub("\\[|\\]|'", "", AA_results$All_genes), ",\\s*")
)

###########################################################
# Detection and location of amplified oncogenes
###########################################################

# number of oncogenes in amplicons
plot_oncogene_counts(AA_results)

# subset for amplicons with amplified oncogenes / amplified genes
amplified_onco <- AA_results[AA_results$Oncogenes != "[]",] # 2139
amplified_gene <- AA_results[AA_results$All_genes != "[]",] # 5864

# amplicon class distribution (amplicons with amplified oncogenes)
plot_amplicon_class(amplified_onco, "amplified_oncogenes")
plot_amplicon_class(amplified_gene, "amplified_all_genes")

###########################################################
# Feature CN across amplicon classes
###########################################################

# helper function to plot CN across amplicon classes
helper_plot_fmcn <- function(amplified, label) {

    # feature median CN AOV and tukey test
    aov.med <- aov(Feature_median_copy_number ~ Classification, data = amplified)
    tukey <- TukeyHSD(aov.med)$Classification |> as.data.frame()
    med_tukey <- tukey[tukey$"p adj" < 0.05,]
    print(med_tukey)

    # feature maximum CN AOV and tukey test
    aov.max <- aov(Feature_maximum_copy_number ~ Classification, data = amplified)
    tukey <- TukeyHSD(aov.max)$Classification |> as.data.frame()
    max_tukey <- tukey[tukey$"p adj" < 0.05,]
    print(max_tukey)

    plot_fmcn_genes(amplified, med_tukey, max_tukey, label)
}

helper_plot_fmcn(amplified_onco, "oncogene")
helper_plot_fmcn(amplified_gene, "all_genes")

###########################################################
# Proportion of ecDNA amplicons carrying amplified oncogenes
###########################################################

ecDNA_amplicons <- AA_results[AA_results$Classification == "ecDNA",]

# plot presence of genes on ecDNA
plot_ecDNA_genes(ecDNA_amplicons, "oncogene")
plot_ecDNA_genes(ecDNA_amplicons, "all")

# ----------------------------------------------------------------------
# Gene Analysis
# ----------------------------------------------------------------------

###########################################################
# Prevalence of amplified oncogenes of ecDNA
###########################################################

# keep only ecDNA amplicons with amplified oncogene
col <- "Oncogenes"
ecDNA_amplicons$amplification <- ifelse(
        ecDNA_amplicons[[col]] == "[]",
        0, 1
    )
amplified_amplicons <- ecDNA_amplicons[ecDNA_amplicons$amplification == 1,]

# get unique oncogenes
oncogenes <- unlist(
  strsplit(gsub("\\[|\\]|'", "", amplified_amplicons$Oncogenes), ",\\s*")
) |> unique()
# logic check
for (gene in oncogenes) {
  if (length(grep(gene, oncogenes)) > 1) print(paste("Need to check", gene))
}

# initate dataframe to store results
toPlot <- data.frame(matrix(nrow=nrow(amplified_amplicons), ncol=length(oncogenes)))
colnames(toPlot) <- c(oncogenes)
rownames(toPlot) <- amplified_amplicons$'Feature ID'

# binarize oncogene detection
for (gene in oncogenes) {
  toPlot[[gene]][grep(gene, amplified_amplicons$Oncogenes)] <- 1
}

# format for plotting
toPlot[is.na(toPlot)] <- 0
toPlot <- t(toPlot) |> as.data.frame()

# plot heatmap
plot_oncogene_heatmap(toPlot)

###########################################################
# Oncogenes amplified in ecDNA vs other classes
###########################################################

# get amplified oncogenes per class
count_bfb <- get_amplified_oncogenes(AA_results, "BFB")
count_cnc <- get_amplified_oncogenes(AA_results, "Complex-non-cyclic")
count_ecd <- get_amplified_oncogenes(AA_results, "ecDNA")
count_lin <- get_amplified_oncogenes(AA_results, "Linear")

# merge data frame
common_genes <- c(
  count_bfb$Gene, count_cnc$Gene, count_ecd$Gene, count_lin$Gene
) |> unique()
toPlot <- data.frame(
  Gene = common_genes,
  BFB = NA,
  CNC = NA,
  ecDNA = NA,
  Linear = NA
)
toPlot$BFB <- count_bfb$BFB[match(toPlot$Gene, count_bfb$Gene)]
toPlot$CNC <- count_cnc$'Complex-non-cyclic'[match(toPlot$Gene, count_cnc$Gene)]
toPlot$ecDNA <- count_ecd$ecDNA[match(toPlot$Gene, count_ecd$Gene)]
toPlot$Linear <- count_lin$Linear[match(toPlot$Gene, count_lin$Gene)]

# format dataframe for plotting
toPlot[is.na(toPlot)] <- 0
toPlot <- toPlot[order(toPlot$ecDNA, toPlot$BFB, toPlot$CNC, toPlot$Linear, decreasing = TRUE),]
rownames(toPlot) <- toPlot$Gene
toPlot$Gene <- NULL
toPlot <- t(toPlot) |> as.data.frame()
toPlot <- toPlot[,-which(colSums(toPlot) == 1)]

# plot heatmap of oncogene detection
plot_class_oncogene_heatmap(toPlot)