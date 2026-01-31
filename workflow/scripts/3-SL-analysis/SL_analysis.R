# load libraries
suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(data.table)
    library(ComplexHeatmap)
    library(RColorBrewer)
})

source("workflow/scripts/utils/palettes.R")

set.seed(101)

###########################################################
# Load in data
###########################################################

# get SL scores
sl_scores <- read.csv("data/results/data/SLscores/BTC.csv")

# get mutations
mut <- fread("workflow/jackie/cbio_mafs/BTC_mutations_extended.txt")
mut <- mut[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Type", "oncogenic")]

# list of oncogenes from OncoKB
oncokb <- fread("data/rawdata/metadata/cancerGeneList.tsv")

###########################################################
# Format dataframe
###########################################################

# create binary matrix
bmat <- mut %>%
    select(-Variant_Type) %>%
    distinct(Tumor_Sample_Barcode, Hugo_Symbol) %>%
    mutate(mutation = 1) %>%
    pivot_wider(
        names_from = Hugo_Symbol,
        values_from = mutation,
        values_fill = 0
    ) %>%
    column_to_rownames("Tumor_Sample_Barcode")

# compute mutation prevalence, keep top 100 prevalent mutations
mut_prevalence <- data.frame(
    Hugo_Symbol = colnames(bmat),
    Count = colSums(bmat)
)
mut_prevalence <- mut_prevalence[order(mut_prevalence$Count, decreasing =  TRUE),]
mut_prevalence <- mut_prevalence[1:100,]

# keep only prevalent mutations
bmat <- bmat[,mut_prevalence$Hugo_Symbol]

# label variant types of sample-gene pairs
variant_anno <- mut %>%
    group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
    summarise(
        Variant_Type = paste(unique(Variant_Type), collapse = " & "),
        .groups = "drop"
    )
# catch duplicates
variant_anno$Variant_Type[variant_anno$Variant_Type == "DEL & SNP"] <- "SNP & DEL"
variant_anno$Variant_Type[variant_anno$Variant_Type == "DNP & SNP"] <- "SNP & DNP"
variant_anno$Variant_Type[variant_anno$Variant_Type == "INS & SNP"] <- "SNP & INS"
variant_anno$Variant_Type[variant_anno$Variant_Type == "DEL & INS"] <- "ING & DEL"


# label bmat with variant type
bmat <- bmat %>%
    rownames_to_column("Tumor_Sample_Barcode") %>%
    pivot_longer(
        cols = -Tumor_Sample_Barcode,
        names_to = "Hugo_Symbol",
        values_to = "mutation"
    ) %>%
    left_join(variant_anno, by = c("Tumor_Sample_Barcode", "Hugo_Symbol")) %>%
    mutate(value = ifelse(mutation == 1, Variant_Type, NA)) %>%
    select(-mutation, -Variant_Type) %>%
    pivot_wider(
        names_from = Hugo_Symbol,
        values_from = value
    ) %>%
    column_to_rownames("Tumor_Sample_Barcode")

# create gene annotations
anno <- data.frame(
    Hugo_Symbol = colnames(bmat),
    Annotation = oncokb$'Gene Type'[match(colnames(bmat), oncokb$'Hugo Symbol')]
)

###########################################################
# Format dataframe
###########################################################

bmat <- t(bmat) |> as.data.frame()

row_ha <- rowAnnotation(
    Variant_Type = variant_anno$Variant_Type[variant_anno$Hugo_Symbol %in% rownames(bmat)],
    col = list(Variant_Type = variant_pal))

# try plot
ht <- Heatmap(
    bmat, name = "Mutations",
    show_column_names = FALSE,
    col = variant_pal,
    na_col = "white",
    rect_gp = gpar(col = "grey80", lwd = 0.5)
)

filename <- paste0("data/results/figures/SL_exploration/heatmap.png")
png(filename, width = 13, height = 15, res = 600, units = "in")
print(draw(ht))
dev.off()


