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

# get SL scores (keep only singleKO)
sl_scores <- read.csv("data/results/data/SLscores/BTC.csv")

# get mutations
mut <- fread("data/procdata/cbio_mafs/BTC_mutations_extended.txt")
mut <- mut[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Type", "oncogenic")]

# list of oncogenes from OncoKB
oncokb <- fread("data/rawdata/metadata/cancerGeneList.tsv")

# list of targeted genes
targets <- read_excel("data/rawdata/metadata/targeted_genes.xlsx", sheet = 1)

###########################################################
# Keep SKOs where mutated gene is a tumour mutations
###########################################################

# keep only SKOs and annotate the mutated gene
sko <- sl_scores[sl_scores$outcome == "singleKO",]
sko$mutated_gene <- ifelse(
    sko$geneA_state == "normal",
    sko$geneB,
    sko$geneA
)
sko$unmutated_gene <- ifelse(
    sko$mutated_gene == sko$geneA,
    sko$geneB,
    sko$geneA
)

# label if mutated gene is a tumour mutation
sko$sample_gene <- paste0(sko$sample, ":", sko$mutated_gene)
mut$sample_gene <- paste0(mut$Tumor_Sample_Barcode, ":", mut$Hugo_Symbol)

sko$tumour_mut <- ifelse(
    sko$sample_gene %in% mut$sample_gene,
    TRUE, FALSE
)

# keep only tumour mutated SKOs
sko_tmut <- sko[sko$tumour_mut == TRUE,]

###########################################################
# Annotate tumour mutated SKOs
###########################################################

# annotate tumour mutations with OncoKB
sko_tmut$mut_type <- oncokb$'Gene Type'[match(sko_tmut$mutated_gene, oncokb$'Hugo Symbol')]

# annotate unmutated gene if targetable
sko_tmut$targetable <- ifelse(sko_tmut$unmutated_gene %in% targets$'Gene Target', TRUE, FALSE)

###########################################################
# Plot count of SKOs
###########################################################

# format dataframe and get counts
count_sko <- sko_tmut %>%
    group_by(sample) %>%
    summarise(
        count_SKO = n(),
        count_targetable = sum(targetable == TRUE)) %>%
    pivot_longer(
        cols = c(count_SKO, count_targetable),
        names_to = "Variable",
        values_to = "Count"
    )
count_sko <- count_sko[order(count_sko$Count, decreasing = TRUE),]
count_sko$sample <- factor(count_sko$sample, levels = rev(unique(count_sko$sample)))

# rename facets
count_sko$Variable[count_sko$Variable == "count_SKO"] <- "All SKOs"
count_sko$Variable[count_sko$Variable == "count_targetable"] <- "Targetable Unmutated Genes"

# plot count of SKOs
ggplot(count_sko, aes(y = sample, x = Count)) +
    geom_bar(stat = "Identity", fill = teal) +
    facet_wrap(~Variable, scales = "free_x") +
    theme_minimal() +
    theme(
        axis.text.y = element_blank()
    ) +
    labs(y = "Patient")


###########################################################
# Top Targetable SKOs
###########################################################

# add in SKO pair name
sko_tmut <- sko_tmut %>%
    mutate(pair = paste(sko_tmut$mutated_gene, sko_tmut$unmutated_gene, sep = " & "))

# get counts of top SKOs
count_tar_sko <- sko_tmut %>%
    filter(targetable == TRUE) %>%
    group_by(pair) %>%
    summarise(Count = n()) %>% 
    unique()
count_tar_sko <- count_tar_sko[order(count_tar_sko$Count, decreasing = TRUE),]
count_tar_sko$pair <- factor(count_tar_sko$pair, levels = rev(count_tar_sko$pair))
count_tar_sko <- count_tar_sko[count_tar_sko$Count > 1,]

# plot count of top targetable SKOs
ggplot(count_tar_sko, aes(y = pair, x = Count)) +
    geom_bar(stat = "Identity", fill = teal) +
    theme_minimal() +
    scale_x_continuous(
        lim = c(0, max(count_tar_sko$Count + 0.5)),
        expand = c(0,0)
    ) +
    labs(y = "Synthetic Lethal Pair")

###########################################################
# Prevalence of Targetable SKOs
###########################################################

# create binary matrix of targetable SKOs
bin_sko <- sko_tmut %>%
    filter(targetable == TRUE) %>%
    filter(pair %in% count_tar_sko$pair) %>%
    select(sample, pair) %>%
    mutate(SKO = 1) %>%
    pivot_wider(
        names_from = pair,
        values_from = SKO,
        values_fill = 0
    ) %>%
    column_to_rownames("sample")

col_ha <- columnAnnotation(
    "Mutated Gene" = sko_tmut$mutated_gene[match(colnames(bin_sko), sko_tmut$pair)],
    "Unmutated Gene" = sko_tmut$unmutated_gene[match(colnames(bin_sko), sko_tmut$pair)],
    "Mutation Type" = sko_tmut$mut_type[match(colnames(bin_sko), sko_tmut$pair)]
)

# plot heatmap
ht <- Heatmap(
    bin_sko, name = "SKO",
    top_annotation = col_ha,
    show_column_names = FALSE,
    show_row_names = FALSE,
    col = rev(binary_pal2),
    na_col = "white",
    rect_gp = gpar(col = "grey80", lwd = 0.5)
)

#filename <- paste0("data/results/figures/SL_exploration/heatmap.png")
#png(filename, width = 13, height = 15, res = 600, units = "in")
print(draw(ht))
#dev.off()

###########################################################
# Prevalence of Targetable Genes
###########################################################

# create binary matrix of targetable genes
bin_tar <- sko_tmut %>%
    mutate(sample_unmut_gene = paste(sko_tmut$sample, sko_tmut$unmutated_gene, sep = ":")) %>%
    filter(targetable == TRUE) %>%
    group_by(sample_unmut_gene) %>%
    summarise(Count = n()) %>%
    mutate(
        sample = sub(":.*", "", sample_unmut_gene),
        gene = sub(".*:", "", sample_unmut_gene)
    ) %>%
    select(-sample_unmut_gene) %>%
    pivot_wider(
        names_from = sample,
        values_from = Count,
        values_fill = 0
    ) %>%
    column_to_rownames("gene")

# get pal
pal <- colorRampPalette(c("#C0CDDD", "#00102F"))(99)
pal <- c("white", pal)

# plot heatmap
ht <- Heatmap(
    bin_tar, name = "No.\nSKOs",
    show_column_names = FALSE,
    col = pal,
    rect_gp = gpar(col = "grey80", lwd = 0.5)
)

draw(ht)
