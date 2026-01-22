# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(dplyr)
})

PROCDATA_DIR <- "data/procdata/"

###########################################################
# Load in data
###########################################################

# load in ecDNA speciifc results
ecDNA_amplicons <- readRDS(paste0(PROCDATA_DIR, "data/ecDNA_amplicons.rds"))

# load in ecDNA counts
ecDNA_counts <- read.csv(paste0(PROCDATA_DIR, "data/ecDNA_counts.csv"))


###########################################################
# Count of ecDNA
###########################################################

toPlot <- ecDNA_counts %>%
    count(cohort, count) %>%
    group_by(cohort) %>%
    mutate(prop = n / sum(n) * 100)
toPlot$count <- factor(toPlot$count, levels = c(max(toPlot$count):0))

ggplot(toPlot, aes(x = cohort, y = prop, fill = factor(count))) +
    geom_col(color = "black") +
    geom_text(
        aes(label = paste0(n, " (", round(prop, 2), "%)")),
        position = position_stack(vjust = 0.5),
        color = "black",
        size = 3
    ) +
    theme_classic() +
    labs(x = "Cohort", y = "Proportion of samples (%)", fill = "ecDNA count\nin sample")
