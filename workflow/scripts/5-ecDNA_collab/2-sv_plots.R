##### quickly check associations

# load libraries
suppressPackageStartupMessages({
    library(data.table)
})

# -- Run on H4H

setwd("/cluster/projects/bhklab/rawdata/SL_MOHCCN/STAR_fusion")
#fusions <- list.files(pattern = "fusions.tsv")
#fusions <- list.files(pattern = "mavis_summary.tab")

res <- data.frame(matrix(nrow=0, ncol=2))
colnames(res) <- c("Sample", "Fusion")
for (file in fusions) {
    fusion_status <- 0
    df <- fread(file, data.table = FALSE)
    colnames(df)[1] <- "gene1"
    f1 <- df[df$gene1 == "MYC" & df$gene2 == "PVT1",]
    f2 <- df[df$gene1 == "PVT1" & df$gene2 == "MYC",]

    if (nrow(f1) > 0 | nrow(f2) > 0) fusion_status <- 1
    res <- rbind(res, data.frame(Sample = file, Fusion = fusion_status))
}

save(res, file = "fusion_res.RData")


# using STAR fusion

fusions <- list.files(pattern = "fusion_predictions.abridged.tsv")

res <- data.frame(matrix(nrow=0, ncol=2))
colnames(res) <- c("Sample", "Fusion")
for (file in fusions) {
    fusion_status <- 0
    df <- fread(file, data.table = FALSE)
    colnames(df)[1] <- "FusionName"
    f1 <- df[df$FusionName == "PVT1--MYC",]
    f2 <- df[df$FusionName == "MYC--PVT1",]

    if (nrow(f1) > 0 | nrow(f2) > 0) fusion_status <- 1
    res <- rbind(res, data.frame(Sample = file, Fusion = fusion_status))
}


###########################################################
# General check
###########################################################

