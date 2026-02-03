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
# Load in data
###########################################################

# get top 50 targetable SKOs
top_SKO_tar <- readRDS("data/results/data/top_SKO_tar.RDS")

# load in CCLE RNAseq
ccle <- fread("data/procdata/psets/CCLE_RSEM_TPM_reheadered.tsv", data.table = FALSE)
rownames(ccle) <- ccle$Hugo_Symbol
ccle$Hugo_Symbol <- NULL

# load in CCLE drug sensitivity
sen <- read.csv("data/procdata/psets/CCLE_sensitivity.csv")
rownames(sen) <- sen$X
sen$X <- NULL

###########################################################
# Check SL status
###########################################################
rna <- ccle
stat_options <- c('1st Qu.', 'Median', 'thirds')
drug <- "Crizotinib"
stat <- stat_options[2]

for (i in 1:nrow(top_SKO_tar)) {

    gm <- top_SKO_tar$mutated_gene[i]
    gu <- top_SKO_tar$unmutated_gene[i]
    
    # subset for SKO pair
    sub <- rna[rownames(rna) %in% c(gm, gu),] |>
        t() |> as.data.frame()
    
    # get summary stats
    stats <- data.frame(
        Mutated = as.numeric(summary(sub[[gm]])),
        Unmutated = as.numeric(summary(sub[[gu]]))
    )
    rownames(stats) <- names(summary(sub[[gm]]))
    stats <- t(stats) |> as.data.frame()
    stats$thirds <- (stats$Max - stats$Min)/3

    # use selected stat to annotate KO
    sub[[gm]] <- ifelse(sub[[gm]] <= stats[[stat]][rownames(stats) == "Mutated"], 0, 1)
    sub[[gu]] <- ifelse(sub[[gu]] <= stats[[stat]][rownames(stats) == "Unmutated"], 0, 1)

    # label SKO
    sub <- sub[complete.cases(sub),]
    sub$status <- NA
    for (i in 1:nrow(sub)) {
        if (sub[[gm]][i] == 0 & sub[[gu]][i] == 0) sub$status[i] <- "DKO"
        if (sub[[gm]][i] == 1 & sub[[gu]][i] == 0) sub$status[i] <- "SKO_targetable"
        if (sub[[gm]][i] == 0 & sub[[gu]][i] == 1) sub$status[i] <- "SKO_other"
        if (sub[[gm]][i] == 1 & sub[[gu]][i] == 1) sub$status[i] <- "noKO"
    }

    # make dataframe of common cells
    common_cells <- intersect(rownames(sub), rownames(sen))
    df <- data.frame(
        sample = common_cells,
        SL = sub$status[match(common_cells, rownames(sub))],
        sen = sen[[drug]][match(common_cells, rownames(sen))]
    )

    ggplot(df, aes(x = SL, y = sen, fill = SL)) +
        geom_boxplot() + geom_jitter(width = 0.2, alpha = 0.5) +
        theme_minimal() +
        labs(y = "AAC", x = "Synthetic Lethal Interaction")



}

###########################################################
# ttest TODO
###########################################################



dko <- sen$Lapatinib[rownames(sen) %in% rownames(sub[sub$status == "DKO",])]
sko <- sen$Lapatinib[rownames(sen) %in% rownames(sub[sub$status != "DKO",])]
t.test(dko, sko, alternative = "two.sided")

