# load libraries
suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(data.table)
    library(ComplexHeatmap)
    library(RColorBrewer)
    library(readxl)
})

source("workflow/scripts/utils/palettes.R")
set.seed(101)

###########################################################
# Load in data
###########################################################

# read in SKO counts
SKO_count <- data.frame(matrix(nrow=0, ncol=0))
for (file in list.files("data/results/data/SL-analysis/SKO_counts", full.names = TRUE)) {
    counts <- read.csv(file)
    SKO_count <- rbind(SKO_count, counts)
}

# read in targetable SKO counts
SKO_tar <- data.frame(matrix(nrow=0, ncol=0))
for (file in list.files("data/results/data/SL-analysis/SKO_targetable", full.names = TRUE)) {
    counts <- read.csv(file)
    SKO_tar <- rbind(SKO_tar, counts)
}

###########################################################
# Format targetable SKO counts
###########################################################

SKO_tar$mutated_gene <- sub(" &.*", "", SKO_tar$pair)
SKO_tar$unmutated_gene <- sub(".*& ", "", SKO_tar$pair)

###########################################################
# Plot SKO counts across cohorts
###########################################################

# boxplot of all SKO counts
filename <- paste0("data/results/figures/SL-analysis/All_Cohorts/All_SKO_counts.png")
png(filename, width = 6, height = 6, res = 600, units = "in")
ggplot(SKO_count[SKO_count$Variable == "All SKOs",], aes(x = Cohort, y = Count)) +
    geom_boxplot(fill = olive) + geom_jitter(width = 0.1, alpha = 0.5) +
    theme_minimal()
dev.off()

# boxplot of all SKO counts
filename <- paste0("data/results/figures/SL-analysis/All_Cohorts/Tar_SKO_counts.png")
png(filename, width = 6, height = 6, res = 600, units = "in")
ggplot(SKO_count[SKO_count$Variable == "Targetable Unmutated Genes",], aes(x = Cohort, y = Count)) +
    geom_boxplot(fill = olive) + geom_jitter(width = 0.1, alpha = 0.5) +
    theme_minimal()
dev.off()

###########################################################
# Plot top targatable SKOs across cohorts
###########################################################

# get top 50 targetable SKOs
top_tar <- SKO_tar %>%
    select(pair, Count) %>%
    group_by(pair) %>%
    mutate(Count_Total = sum(Count)) %>%
    select(-Count) %>%
    unique()
top_tar <- top_tar[order(top_tar$Count_Total, decreasing = TRUE),]
top_tar <- top_tar$pair[1:50]

# keep top 50 targetable SKOs
top_SKO_tar <- SKO_tar[SKO_tar$pair %in% top_tar,]
top_SKO_tar$pair <- factor(top_SKO_tar$pair, levels = rev(top_tar))

# plot top 50 targetable SKOs
filename <- paste0("data/results/figures/SL-analysis/All_Cohorts/top50_tar_SKOs.png")
png(filename, width = 6, height = 6, res = 600, units = "in")
ggplot(top_SKO_tar, aes(x = Count, y = pair, fill = Cohort)) +
    geom_bar(stat = "identity") +
    scale_x_continuous(lim = c(0, 105), expand = c(0,0)) +
    scale_fill_manual(values = cohort_pal) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 6)) +
    labs(y = "Synthetic Lethal Pair")
dev.off()