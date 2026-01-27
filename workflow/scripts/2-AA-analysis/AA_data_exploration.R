# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggridges)
    library(dplyr)
    library(reshape2)
    library(ComplexHeatmap)
    library(RColorBrewer)
})

PROCDATA_DIR <- "data/procdata/"
RESULTS_DIR <- "data/results/"

source("workflow/scripts/utils/plots.R")
source("workflow/scripts/utils/palettes.R")

###########################################################
# Load in data
###########################################################

# load in AA results
AA_results <- readRDS(paste0(PROCDATA_DIR, "AA_output/AA_results.rds"))
AA_results$cohort <- sub("_.*", "", AA_results$'Sample name')
AA_results <- AA_results[!is.na(AA_results$'AA amplicon number'),]

# load in ecDNA speciifc results
ecDNA_amplicons <- readRDS(paste0(PROCDATA_DIR, "AA_output/ecDNA_amplicons.rds"))
ecDNA_amplicons$cohort <- sub("_.*", "", ecDNA_amplicons$'Sample name')

# load in ecDNA counts
ecDNA_counts <- read.csv(paste0(PROCDATA_DIR, "AA_output/ecDNA_counts.csv"))

###########################################################
# Number of samples (w amplicon) per cohort
###########################################################

toPlot <- table(unique(AA_results[,1:2])$cohort) |> as.data.frame()


###########################################################
# Number of amplicons in each sample
###########################################################

# get number of amplicons per sample
num_amplicons <- data.frame(table(AA_results$'Sample name'))
num_amplicons$cohort <- sub("_.*", "", num_amplicons$Var1)

# plot number of amplicons
plot_num_amplicons(num_amplicons)
plot_num_amplicons2(num_amplicons)

###########################################################
# Number of genes in each sample
###########################################################

# get number of amplified oncogenes
AA_results$n_genes <- lengths(
    strsplit(gsub("\\[|\\]|'", "", AA_results$Oncogenes), ",\\s*")
)

# get number of amplified genes
AA_results$n_genes_all <- lengths(
    strsplit(gsub("\\[|\\]|'", "", AA_results$'All genes'), ",\\s*")
)

plot_num_genes(AA_results, "oncogene")
plot_num_genes(AA_results, "all")

###########################################################
# Count of ecDNA per sample
###########################################################

# count of amplicons on ecDNA
plot_ecDNA_counts(ecDNA_counts)

# binary presence of amplicons on ecDNA
plot_ecDNA_counts(ecDNA_counts, type = "binary")

###########################################################
# Detection and location of amplified oncogenes
###########################################################

# number of oncogenes in amplicons
plot_oncogene_counts(AA_results)

# subset for amplicons with amplified oncogenes
amplified <- AA_results[AA_results$Oncogenes != "[]",]

# amplicon class distribution (amplicons with amplified oncogenes)
plot_amplicon_class(amplified)

###########################################################
# Proportion of ecDNA amplicons carrying amplified oncogenes
###########################################################

# plot presence of genes on ecDNA
plot_ecDNA_genes(ecDNA_amplicons, "oncogene")
plot_ecDNA_genes(ecDNA_amplicons, "all")

###########################################################
# Prevalence of amplified oncogenes of ecDNA
###########################################################

# keep only ecDNA amplicons with amplified oncogene
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
