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

###########################################################
# Number of samples (w amplicon) per cohort
###########################################################

plot_num_patients(AA_results)
length(unique(AA_results$'Sample name'[AA_results$Centre == "PM2C"]))
length(unique(AA_results$'Sample name'[AA_results$Centre == "BCCA"]))

###########################################################
# Number of amplicons in each sample
###########################################################

# get number of amplicons per sample
num_amplicons <- data.frame(table(AA_results$'Sample name'))
num_amplicons$cohort <- AA_results$cohort[match(num_amplicons$Var1, AA_results$'Sample name')]
num_amplicons$Centre <- AA_results$Centre[match(num_amplicons$Var1, AA_results$'Sample name')]

# plot number of amplicons
plot_num_amplicons(num_amplicons)

###########################################################
# Distribution of amplicon class among available amplicons
###########################################################

plot_amplicon_class(AA_results, "all_amplicons")

###########################################################
# Distribution of feature copy  number per amplicon class
###########################################################

plot_fmcn(AA_results, "Feature median copy number")
plot_fmcn(AA_results, "Feature maximum copy number")

###########################################################
# Count of ecDNA per patient
###########################################################

ecDNA <- AA_results[AA_results$Classification == "ecDNA",]
ecDNA_counts <- data.frame(table(ecDNA$'Sample name'))
colnames(ecDNA_counts) <- c("sample", "count")
ecDNA_counts$cohort <- ecDNA$cohort[match(ecDNA_counts$sample, ecDNA$'Sample name')]

# count of amplicons on ecDNA
plot_ecDNA_prop(ecDNA_counts)

# binary presence of amplicons on ecDNA (including no ecDNA)
no_ecDNA <- AA_results[-which(AA_results$'Sample name' %in% ecDNA_counts$sample),c(2,1)]
no_ecDNA <- data.frame(sample = no_ecDNA$'Sample name', count = 0, cohort = no_ecDNA$cohort)
ecDNA_counts <- rbind(ecDNA_counts, no_ecDNA)
plot_ecDNA_prop(ecDNA_counts, type = "binary")