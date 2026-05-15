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

source("workflow/scripts/2-AA-analysis/utils/get_tables.R")
source("workflow/scripts/utils/cohorts.R")
source("workflow/scripts/utils/palettes.R")

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
# Check duplicates
###########################################################

dups <- all_samples[-which(all_samples$duplicated == "not_duplicated"),]
toPlot <- table(unique(dups[,c(3:4)])$cohort) |> as.data.frame()
toPlot$Centre <- ifelse(toPlot$Var1 %in% PM2C_cohorts, "PM2C", "BCCA")
toPlot$Centre <- factor(toPlot$Centre, levels = names(centre_pal))
toPlot$Var1 <- factor(toPlot$Var1, levels = c(PM2C_cohorts, BCCA_cohorts))

p <- ggplot(toPlot, aes(x = Var1, y = Freq, fill = Centre)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(
            aes(label = Freq),
            vjust = -0.3,
            color = "black",
            size = 3
        ) +
  scale_fill_manual(values = centre_pal) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Count", x = "Cohort")
filename <- "data/results/figures/AA_exploration/processing/duplicates.png"
ggsave(filename, p)

###########################################################
# Remove duplicates
###########################################################

all_samples <- all_samples[-which(all_samples$duplicated == "duplicated_removed"),]

###########################################################
# Number with no amplicons
###########################################################

no_amps <- all_samples[all_samples$result_table == "missing",] #805
toPlot <- table(unique(no_amps[,c(3:4)])$cohort) |> as.data.frame()
toPlot$Centre <- ifelse(toPlot$Var1 %in% PM2C_cohorts, "PM2C", "BCCA")
toPlot$Centre <- factor(toPlot$Centre, levels = names(centre_pal))
toPlot$Var1 <- factor(toPlot$Var1, levels = c(PM2C_cohorts, BCCA_cohorts))

# regular bar plot
p <- ggplot(toPlot, aes(x = Var1, y = Freq, fill = Centre)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(
            aes(label = Freq),
            vjust = -0.3,
            color = "black",
            size = 3
        ) +
  scale_fill_manual(values = centre_pal) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Count", x = "Cohort")
filename <- "data/results/figures/AA_exploration/processing/no_amps.png"
ggsave(filename, p)

# group by duplication status
toPlot <- unique(no_amps[,c(3:4,9)])
toPlot <- table(toPlot$cohort, toPlot$duplicated) |> as.data.frame()
toPlot$Centre <- ifelse(toPlot$Var1 %in% PM2C_cohorts, "PM2C", "BCCA")
toPlot$Centre <- factor(toPlot$Centre, levels = names(centre_pal))
toPlot$Var1 <- factor(toPlot$Var1, levels = c(PM2C_cohorts, BCCA_cohorts))

p <- ggplot(toPlot, aes(x = Var1, y = Freq, fill = Var2, group = Var2)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge2()) +
  geom_text(
            aes(label = Freq),
            vjust = -0.3,
            color = "black",
            size = 3,
            position = position_dodge2(width = 0.9)
        ) +
  scale_fill_manual("Duplicated", values = binary_pal3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Count", x = "Cohort")
filename <- "data/results/figures/AA_exploration/processing/no_amps_dups.png"
ggsave(filename, p)

###########################################################
# Remove samples with missing output files
###########################################################

all_samples <- all_samples[-which(all_samples$result_table == "missing"),]

###########################################################
# Identify sample mismatches
###########################################################

# 4 mismatched sample ids
table(AA_results$'Sample name' %in% all_samples$sample)
table(all_samples$sample %in% AA_results$'Sample name')

AA_results$'Sample name'[-which(AA_results$'Sample name' %in% all_samples$sample)]
all_samples$sample[-which(all_samples$sample %in% AA_results$'Sample name')]

###########################################################
# Number with NA amplicons
###########################################################

na_amp <- AA_results[is.na(AA_results$'AA amplicon number'),]

toPlot <- table(unique(na_amp[,c(1:2)])$cohort) |> as.data.frame()
toPlot$Centre <- ifelse(toPlot$Var1 %in% PM2C_cohorts, "PM2C", "BCCA")
toPlot$Centre <- factor(toPlot$Centre, levels = names(centre_pal))
toPlot$Var1 <- factor(toPlot$Var1, levels = c(PM2C_cohorts, BCCA_cohorts))

p <- ggplot(toPlot, aes(x = Var1, y = Freq, fill = Centre)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(
            aes(label = Freq),
            vjust = -0.3,
            color = "black",
            size = 3
        ) +
  scale_fill_manual(values = centre_pal) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Count", x = "Cohort")
filename <- "data/results/figures/AA_exploration/processing/na_amp_detected.png"
ggsave(filename, p)