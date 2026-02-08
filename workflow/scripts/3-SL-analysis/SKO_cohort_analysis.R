# load libraries
suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(data.table)
    library(ComplexHeatmap)
    library(RColorBrewer)
    library(readxl)
    library(reshape2)
    library(ggbreak)
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
toPlot$variable <- factor(toPlot$variable, levels = rev(colnames(cohort_count)[2:4]))
toPlot$Cohort <- factor(toPlot$Cohort, levels = cohort_count$Cohort[order(cohort_count$All)])

filename <- paste0("data/results/figures/SL-analysis/All_Cohorts/Cohort_Counts.png")
png(filename, width = 4, height = 3.5, res = 600, units = "in")
ggplot(toPlot, aes(y = Cohort, x = value, fill = variable)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "black", linewidth = 0.2) +
    scale_fill_manual("", values = three_pal) +
    theme_minimal() +
    theme(legend.position = "none", panel.border = element_rect()) +
    labs(x = "Number of Patients")
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

###########################################################
# Counta of targetable SKOs across mutations
###########################################################

# count of top mutated genes
SKO_tar$mutated_gene <- factor(SKO_tar$mutated_gene, levels = unique(SKO_tar$mutated_gene))
count_mut <- SKO_tar %>%
    group_by(mutated_gene) %>%
    summarise(Count = sum(Count))
top_mut <- count_mut[order(count_mut$Count, decreasing = TRUE),]$mutated_gene[1:20]

# count of top unmutated genes
count_unmut <- SKO_tar %>%
    group_by(unmutated_gene) %>%
    summarise(Count = sum(Count))
top_unmut <- count_unmut[order(count_unmut$Count, decreasing = TRUE),]$unmutated_gene

# count of SKOs in patinets across mutated genes
count_SKO <- SKO_tar %>%
    group_by(mutated_gene, unmutated_gene) %>%
    summarise(Count = sum(Count))
count_SKO <- count_SKO[count_SKO$mutated_gene %in% top_mut,]
count_SKO$mutated_gene <- factor(count_SKO$mutated_gene, levels = rev(top_mut))
count_SKO$unmutated_gene <- factor(count_SKO$unmutated_gene, levels = top_unmut)


# plot prevalence by mutated genes (normal axis)
filename <- paste0("data/results/figures/SL-analysis/All_Cohorts/top_mutated_ynormal.png")
png(filename, width = 4.5, height = 4, res = 600, units = "in")
ggplot(count_SKO, aes(x = log(Count), y = mutated_gene, fill = unmutated_gene)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
    scale_fill_manual("Unmutated Gene", values = palette_20) +
    theme_minimal() +
    theme(legend.position = "none", panel.border = element_rect()) +
    labs(y = "Mutated Gene", x = "Log-Normalized Count of Targetable SKOs")
dev.off()

# plot prevalence by mutated genes (break axis)
filename <- paste0("data/results/figures/SL-analysis/All_Cohorts/top_mutated_ybreak.png")
png(filename, width = 6, height = 4, res = 600, units = "in")
ggplot(count_SKO, aes(x = Count, y = mutated_gene, fill = unmutated_gene)) +
    geom_bar(stat = "identity", color = "gray", linewidth = 0.2) +
    scale_x_break(c(75, 550), scale = 0.2, ticklabels = c(550, 600)) +
    scale_x_break(c(620, 900), scale = 0.2, ticklabels = c(900, 950)) +
    scale_fill_manual(values = palette_20) +
    theme_minimal() +
    theme(legend.position = "none", panel.border = element_rect()) +
    labs(y = "Mutated Gene", x = "Count of Targetable SKOs in Patients")
dev.off()

###########################################################
# Plot pie chart
###########################################################

toPlot <- data.frame(
    Label = c("Total", "SKO"),
    Prop = c(100-76.2, 76.2)
)

ggplot(toPlot, aes(x="", y=Prop, fill=Label)) +
  geom_bar(stat="identity", width=1, color="black", linewidth=1.25) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = c("#812a2aff", "#ddd4d4")) +
  theme_void()
