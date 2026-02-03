# load libraries
suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(data.table)
})

source("workflow/scripts/utils/palettes.R")

###########################################################
# Get number of patients per cohort
###########################################################

cohorts <- c(
    "BTC",
    "COMPASS",
    "GOTIT",
    "INSPIRE",
    "IOALINES",
    "IRIS",
    "MOHCCNO",
    "OCTANE"
)
count <- data.frame(Cohort = cohorts, Count = NA)

for (i in 1:nrow(count)) {
    cohort = count$Cohort[i]
    bn <- paste0("data/results/data/binarymat/", cohort, ".csv")
    mat <- fread(bn)
    count$Count[i] <- length(colnames(mat))-1
}

# save counts
saveRDS(count, file = "data/results/data/cohort_counts.RDS")
