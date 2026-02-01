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

# get mutations
mut <- fread("data/procdata/cbio_mafs/BTC_mutations_extended.txt")
mut <- mut[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Type", "oncogenic")]

# list of oncogenes from OncoKB
oncokb <- fread("data/rawdata/metadata/cancerGeneList.tsv")

###########################################################
# Count of gene mutations
###########################################################

# count per gene-variant type
variant_count <- mut %>%
    group_by(Hugo_Symbol, Variant_Type) %>%
    summarise(Count = n())

# count per gene
gene_count <- mut %>%
    group_by(Hugo_Symbol) %>%
    summarise(Count = n())

# keep top 50 mutated genes
gene_count <- gene_count[order(gene_count$Count, decreasing = TRUE),]
gene_count <- gene_count[1:50,]
variant_count <- variant_count[variant_count$Hugo_Symbol %in% gene_count$Hugo_Symbol,]

# formatting for plotting
variant_count$Hugo_Symbol <- factor(variant_count$Hugo_Symbol, levels = rev(gene_count$Hugo_Symbol))
variant_count$Variant_Type <- factor(variant_count$Variant_Type, levels = names(variant_pal))

###########################################################
# Plot top mutations
###########################################################

# create gene annotations
anno <- data.frame(
    Hugo_Symbol = variant_count$Hugo_Symbol,
    Annotation = oncokb$'Gene Type'[match(variant_count$Hugo_Symbol, oncokb$'Hugo Symbol')]
)

p <- ggplot(variant_count, aes(y = Hugo_Symbol, x = Count, fill = Variant_Type)) + 
    geom_bar(stat = "Identity") +
    scale_fill_manual(values = variant_pal) +
    theme_minimal() +
    scale_x_continuous(lim = c(0, max(variant_count$Count + 2)), expand = c(0,0))
