# load libraries
suppressPackageStartupMessages({
    library(PharmacoGx)
    library(tidyverse)
    library(ggplot2)
    library(data.table)
})

source("workflow/scripts/utils/palettes.R")
set.seed(101)

###########################################################
# Load in data (TODO:: MOVE THIS TO 0-)
###########################################################

# get drug sensitivity
sen <- PharmacoGx::summarizeSensitivityProfiles(
    ccle,
    sensitivity.measure = "aac_recomputed",
    fill.missing = FALSE,
    summary.stat = "median"
) |> t() |> as.data.frame()



###########################################################
# Check SL status
###########################################################

# subset for SKO pair
#sub <- rna[c(635, 8125),]
sub <- rna[c(8866, 8125),]
#sub <- rna[c(8866, 6783),]
rownames(sub) <- c("ROS1", "TP53")
sub <- t(sub) |> as.data.frame()

stats <- data.frame(
    ROS1 = as.numeric(summary(sub$ROS1)),
    TP53 = as.numeric(summary(sub$TP53))
)
rownames(stats) <- names(summary(sub$ROS1))
stats <- t(stats) |> as.data.frame()
stats$thirds <- (stats$Max - stats$Min)/3

# try q1 first
sub$ROS1 <- ifelse(sub$ROS1 <= stats$'1st Qu.'[rownames(stats) == "ROS1"], 0, 1)
sub$TP53 <- ifelse(sub$ROS1 <= stats$'1st Qu.'[rownames(stats) == "TP53"], 0, 1)

# label SKO
sub <- sub[complete.cases(sub),]
sub$status <- NA
for (i in 1:nrow(sub)) {
    if (sub$ROS1[i] == 0 & sub$TP53[i] == 0) sub$status[i] <- "DKO"
    if (sub$ROS1[i] == 1 & sub$TP53[i] == 0) sub$status[i] <- "SKO1"
    if (sub$ROS1[i] == 0 & sub$TP53[i] == 1) sub$status[i] <- "SKO2"
    if (sub$ROS1[i] == 1 & sub$TP53[i] == 1) sub$status[i] <- "noKO"
}


###########################################################
# Plot Crizotinib response
###########################################################

common_cells <- intersect(rownames(sub), rownames(sen))



df <- data.frame(
    sample = common_cells,
    SL = sub$status[match(common_cells, rownames(sub))],
    sen = sen$Lapatinib[match(common_cells, rownames(sen))]
)


ggplot(df, aes(x = SL, y = sen, fill = SL)) +
    geom_boxplot() + geom_jitter(width = 0.1, alpha = 0.5) +
    theme_minimal() +
    scale_fill_manual(values = c(olive, teal)) +
    labs(y = "AAC", x = "Synthetic Lethal Interaction")


dko <- sen$Lapatinib[rownames(sen) %in% rownames(sub[sub$status == "DKO",])]
sko <- sen$Lapatinib[rownames(sen) %in% rownames(sub[sub$status != "DKO",])]
t.test(dko, sko, alternative = "two.sided")

Gene: KRAS (ENSG00000133703) - Summary
