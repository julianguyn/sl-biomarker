##### quickly check associations

# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(dplyr)
    library(purrr)
    library(broom)
    library(ggrepel)
})

PROCDATA_DIR <- "data/procdata/"
RESULTS_DIR <- "data/results/"

source("workflow/scripts/utils/palettes.R")
source("workflow/scripts/utils/cohorts.R")
source("workflow/scripts/2-AA-analysis/utils/utils.R")

###########################################################
# Load in data
###########################################################

# load in metadata
PM_meta <- read.csv("metadata/2026-05-29_PM_meta.csv")
BC_meta <- read.csv("metadata/2026-05-29_BC_meta.csv")

# load in RNA
PM_rna <- readRDS("data/procdata/RNA_TPM/PM_RNA.rds")

# load in mutations
PM_mut <- read.table("data/procdata/mutations/PM_mutations.tsv")
rownames(PM_mut) <- gsub("\\.", "-", rownames(PM_mut))

# load in AA results
AA_results <- readRDS(paste0(PROCDATA_DIR, "AA_output/AA_results.rds"))

# load in all samples
all_samples <- readRDS(paste0(PROCDATA_DIR, "AA_output/sample_df.rds"))

# process cohort names (from utils/get_tables.R)
AA_results <- process_cohort_names(AA_results, res = TRUE)
all_samples <- process_cohort_names(all_samples)

###########################################################
# Preprocessing
###########################################################

# remove duplicates
all_samples <- all_samples[-which(all_samples$duplicated == "duplicated_removed"),]

# remove missing amplicons
all_samples <- all_samples[-which(all_samples$result_table == "missing"),]

# remove NA amplicons
AA_results <- AA_results[!is.na(AA_results$'AA amplicon number'),]
all_samples <- all_samples[all_samples$sample %in% AA_results$'Sample name',]

# add patient ID column
AA_results$Patient_ID <- sub("^([A-Za-z0-9]+_[0-9]+)_.*$", "\\1", AA_results$'Sample name')

# get ecDNA
ecDNA <- AA_results[
    AA_results$Classification == "ecDNA",
    c("Centre", "cohort", "Sample name", "Patient_ID", "Feature ID", "Oncogenes", "Feature median copy number", "Feature maximum copy number")
] |> as.data.frame()

###########################################################
# Make ecDNA - mutation meta for all samples
###########################################################

ecDNA_mut <- data.frame(sample = all_samples$sample)
ecDNA_mut$ecDNA <- ifelse(ecDNA_mut$sample %in% ecDNA$'Sample name', "ecDNA+", "ecDNA-")
ecDNA_mut$match_id <- sub("_WG", "", ecDNA_mut$sample)

# 1213 samples in toPlot

###########################################################
# Plot TP53 - ecDNA status proportion plot
###########################################################

# keep samples with TP53 mutation data
toPlot <- ecDNA_mut[ecDNA_mut$match_id %in% rownames(PM_mut),]
toPlot$TP53 <- PM_mut$TP53[match(toPlot$match_id, rownames(PM_mut))]
toPlot$TP53 <- ifelse(toPlot$TP53 > 0, "Mut", "Wt")

# 968 samples

# get counts
totals <- toPlot %>%
  count(ecDNA, name = "total")

# get proportions
toPlot <- toPlot %>%
    count(ecDNA, TP53) %>%
    group_by(ecDNA) %>%
    mutate(prop = n / sum(n))
toPlot$TP53 <- factor(toPlot$TP53, levels = c("Wt", "Mut"))

# plot
p <- ggplot(toPlot, aes(x = ecDNA, y = prop, fill = TP53)) +
    geom_col(position = "fill", width = 0.6, color = "black") +
    geom_text(
        data = totals,
        aes(x = ecDNA, y = 1.05, label = paste0("n = ", total)),
        inherit.aes = FALSE,
        size = 4
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual("TP53 status", values = c("Mut" = "#469D77", "Wt" = "gray")) +
    labs(x = "", y = "% Tumor Samples", fill = "TP53") +
    theme_classic()
filename <- "data/results/figures/ecDNA_collab/proportion_plots.png"
ggsave(filename, p, w=4, h=3)

###########################################################
# Plot common OR volcano plot
###########################################################

# add ecDNA to PM_mut
genes <- colnames(PM_mut)
PM_mut$ecDNA <- ecDNA_mut$ecDNA[match(rownames(PM_mut), ecDNA_mut$match_id)]
PM_mut <- PM_mut %>%
  mutate(across(-ecDNA, ~ if_else(. > 0, "Mut", "Wt")))


test_genes <- genes[1:2000]
# loop Fisher's test
res <- map_dfr(test_genes, function(gene) {
  df <- PM_mut %>%
    dplyr::filter(!is.na(.data[[gene]]), .data[[gene]] %in% c("Mut", "Wt")) %>%
    dplyr::mutate(status = factor(.data[[gene]], levels = c("Mut", "Wt")))
  
  tab <- table(df$ecDNA, df$status)
  
  # skip genes that don't form 2x2 after filtering
  if (!all(dim(tab) == c(2, 2))) {
    return(tibble(gene = gene, OR = NA, p.value = NA, ci_low = NA, ci_high = NA))
  }
  ft <- fisher.test(tab)
  tibble(
    gene = gene,
    OR = ft$estimate,
    p.value = ft$p.value,
    ci_low = ft$conf.int[1],
    ci_high = ft$conf.int[2]
  )
})

#res <- readRDS("data/procdata/mutations/OR.RDS")

toPlot <- res %>%
    mutate(
    log2OR = -log2(OR),
    padj = p.adjust(p.value, method = "BH"),
    sig = case_when(
      padj < 0.05 & log2OR > 0 ~ "Enriched",
      padj < 0.05 & log2OR < 0 ~ "Depleted",
      TRUE ~ "NS"
    ))

# drop the inf ORs
toPlot <- toPlot %>%
  mutate(
    OR = case_when(
      is.infinite(OR) ~ NA_real_,
      TRUE ~ OR
    ),
    log2OR = ifelse(is.infinite(log2OR), NA_real_, log2OR)
  )

# plot volcano plot
quick_pal <- c("Depleted" = "#8ABEE1", "Enriched" = "#469D77", NS = "gray")
p <- ggplot(toPlot, aes(x = log2OR, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.6) +
    geom_text_repel(
        data = toPlot[toPlot$gene == "TP53",],
        aes(label = gene),
        size = 3.5,
        show.legend = FALSE
    ) +
    scale_color_manual("", values = quick_pal, labels = c("ecDNA depleted", "ecDNA enriched", "Not significant")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_classic() +
    labs(x = "log2(Common Odds Ratio)", y = "-log10(FDR)")
filename <- "data/results/figures/ecDNA_collab/volcano_plot.png"
ggsave(filename, p, w=5, h=3)

