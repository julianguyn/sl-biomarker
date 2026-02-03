# load libraries
suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(data.table)
    library(ComplexHeatmap)
    library(RColorBrewer)
    library(readxl)
    library(reshape2)
})

source("workflow/scripts/utils/palettes.R")
set.seed(101)

###########################################################
# Load in data
###########################################################

# load in cohort counts
cohort_count <- readRDS("data/results/data/cohort_counts.RDS")

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
# Number of patients with SKOs
###########################################################

# get number of SKOs
count_SKO <- SKO_count %>%
    filter(Count > 0) %>%
    select(-Count) %>%
    unique() %>%
    group_by(Cohort, Variable) %>%
    summarise(Count_SKO = n()) %>%
    pivot_wider(
        names_from = Variable,
        values_from = Count_SKO
    )

# get number of targetable SKOs
cohort_count <- cbind(cohort_count, count_SKO) |> as.data.frame()
cohort_count[,3] <- NULL

# get proportion
cohort_count$Prop_SKOs <- cohort_count$'All SKOs' / cohort_count$Count * 100
cohort_count$Prop_tar <- cohort_count$'Targetable Unmutated Genes' / cohort_count$Count * 100

###########################################################
# Plot cohort counts
###########################################################

colnames(cohort_count)[2:4] <- c("All", "All SKOs", "Targetable SKOs")
toPlot <- reshape2::melt(cohort_count)
toPlot <- toPlot[-which(toPlot$variable %in% c("Prop_SKOs", "Prop_tar")),]

filename <- paste0("data/results/figures/SL-analysis/All_Cohorts/Cohort_Counts.png")
png(filename, width = 7, height = 5, res = 600, units = "in")
ggplot(toPlot, aes(x = Cohort, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual("", values = c(teal, olive, ash)) +
    theme_minimal() +
    labs(y = "Number of Patients")
dev.off()

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

# save top targetable SKOs
saveRDS(top_SKO_tar, file = "data/results/data/top_SKO_tar.RDS")


