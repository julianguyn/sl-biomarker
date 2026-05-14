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
AA_results$cohort <- sub("_.*", "", AA_results$'Sample name')
AA_results <- AA_results[!is.na(AA_results$'AA amplicon number'),]

# load in ecDNA speciifc results
ecDNA_amplicons <- readRDS(paste0(PROCDATA_DIR, "AA_output/ecDNA_amplicons.rds"))
ecDNA_amplicons$cohort <- sub("_.*", "", ecDNA_amplicons$'Sample name')

# load in ecDNA counts
ecDNA_counts <- read.csv(paste0(PROCDATA_DIR, "AA_output/ecDNA_counts.csv"))

###########################################################
# Process cohort names
###########################################################

# process AA_results
PM2C_AA_results <- AA_results[AA_results$cohort %in% PM2C_cohorts,]
BCCA_AA_results <- AA_results[-which(AA_results$cohort %in% PM2C_cohorts),]

for (cohort in BCCA_cohorts) {
  BCCA_AA_results$cohort[grep(cohort, BCCA_AA_results$cohort)] <- cohort
}

PM2C_AA_results$Centre <- "PM2C"
BCCA_AA_results$Centre <- "BCCA"

AA_results <- rbind(PM2C_AA_results, BCCA_AA_results)

# process ecDNA counts
ecDNA_counts$cohort <- AA_results$cohort[match(ecDNA_counts$sample, AA_results$'Sample name')]
ecDNA_counts$Centre <- AA_results$Centre[match(ecDNA_counts$sample, AA_results$'Sample name')]

###########################################################
# Number of samples (w amplicon) per cohort
###########################################################

plot_num_patients(AA_results)

###########################################################
# Number of amplicons in each sample
###########################################################

# get number of amplicons per sample
num_amplicons <- data.frame(table(AA_results$'Sample name'))
num_amplicons$cohort <- AA_results$cohort[match(num_amplicons$Var1, AA_results$'Sample name')]
num_amplicons$Centre <- AA_results$Centre[match(num_amplicons$Var1, AA_results$'Sample name')]

# plot number of amplicons
plot_num_amplicons(num_amplicons) # without outliers, highest is 62 (GOLDEN)
#plot_num_amplicons2(num_amplicons)

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

#plot_num_genes(AA_results, "oncogene")
#plot_num_genes(AA_results, "all")

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