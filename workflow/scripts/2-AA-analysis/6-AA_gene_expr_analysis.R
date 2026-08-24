# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

PROCDATA_DIR <- "data/procdata/"
RESULTS_DIR <- "data/results/"

source("workflow/scripts/utils/palettes.R")

###########################################################
# Load in data
###########################################################

oncogene_df <- readRDS(paste0(PROCDATA_DIR, "AA_exploration/oncogene_df.rds"))

# load in RNA
PM_rna <- readRDS("data/procdata/RNA_TPM/PM_RNA.rds")

# load in metadata
PM_meta <- read.csv("metadata/2026-05-29_PM_meta.csv")

# order oncogenes to investigate
count_oncogene <- data.frame(table(oncogene_df$oncogene))
count_oncogene <- count_oncogene[order(count_oncogene$Freq, decreasing = TRUE),]


###########################################################
# ecDNA amplification and gene expression
###########################################################

## TODO: turn into function

gene <- "MYC"
class <- "ecDNA"

gene_exp <- rna_df[rownames(rna_df) == gene,]
amps <- oncogene_df[oncogene_df$oncogene == gene & oncogene_df$Classification == class,]

#' To do
#' Two analyses: 
#' 1) if +ecDNA -> higher expression than -ecDNA
#' 2) if higher CN in ecDNA -> higher expression than lower CN

### 1) preprocessing (will change based on RNA expression matrices)

colnames(gene_exp) <- amps$Sample_name[1:31]

dummy_data <- data.frame(
    MYC = sample(1:5, length(LETTERS), replace = TRUE) * 0.2
)
rownames(dummy_data) <- paste0("patient_", LETTERS)
dummy_data <- t(dummy_data)

gene_exp <- t(cbind(gene_exp, dummy_data)) |> as.data.frame()

gene_exp$label <- c(rep("ecDNA+", 31), rep("ecDNA-", 26))

### 2) ecDNA+ vs ecDNA-

ggplot(gene_exp, aes(x = "MYC", y = .data[[gene]], fill = label)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(shape = 21, alpha = 0.7) +
    scale_fill_manual(values = rev(binary_pal3)) +
    theme_minimal()
