# load libraries
suppressPackageStartupMessages({
    library(PharmacoGx)
    library(tidyverse)
    library(ggplot2)
    library(data.table)
    library(readxl)
})

source("workflow/scripts/utils/palettes.R")
set.seed(101)

###########################################################
# Load in data
###########################################################

# get top 50 targetable SKOs
top_SKO_tar <- readRDS("data/results/data/top_SKO_tar.RDS")

# read in targetable genes
targets <- read_excel("data/rawdata/metadata/targeted_genes.xlsx", sheet = 1)
colnames(targets) <- c("Gene", "Drug", "Brand", "Cancer")
targets$Drug <- str_to_sentence(targets$Drug)

# load in drug sensitivities
load("data/procdata/psets/sensitivity_data.RData")

###########################################################
# TODO:: Turn this into function:
###########################################################

# args:
filepath <- "data/procdata/psets/CCLE_RSEM_TPM_reheadered.tsv"
stat_options <- c('1st Qu.', 'Median', 'thirds')
drug <- "Crizotinib"
stat <- stat_options[1]
sen <- ctrp_sen
pset <- "CTRP2"

# load in rna
rna <- fread(filepath, data.table = FALSE)
rownames(rna) <- rna$Hugo_Symbol
rna$Hugo_Symbol <- NULL


for (i in 1:nrow(top_SKO_tar)) {

    # get SKO pair
    gm <- top_SKO_tar$mutated_gene[i]
    gu <- top_SKO_tar$unmutated_gene[i]

    # get targeted therapies
    gu_targ <- targets$Drug[targets$Gene == gu]
    sen_targ <- sen[rownames(sen) %in% gu_targ,]
    if (nrow(sen_targ) == 0) next
    sen_targ <- t(sen_targ) |> as.data.frame()
    
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
    sub$status <- factor(sub$status, levels = names(status_pal))

    # make dataframe of common cells
    common_cells <- intersect(rownames(sub), rownames(sen_targ))
    df <- data.frame(
        sample = common_cells,
        SL = sub$status[match(common_cells, rownames(sub))]
    )

    # iterate for every available drug
    for (drug in rownames(sen_targ)) {
        df$sensitivity <- sen_targ[[drug]][match(common_cells, rownames(sen_targ))]

        # plot boxplots
        p <- ggplot(df, aes(x = SL, y = sensitivity, fill = SL)) +
        geom_boxplot() + geom_jitter(width = 0.2, alpha = 0.5) +
        scale_fill_manual("Synthetic Lethal\nInteraction", values = SL_pal) +
        theme_minimal() +
        theme(panel.border = element_rect()) +
        labs(y = paste(drug, "Response (AAC)"), x = "Synthetic Lethal Interaction") +
        ggtitle(paste0("Dataset: ", pset, ", SKO: ", gm, " & ", gu))
        
        filename <- paste0("data/results/figures/CCL-top-SKOs/", pset, "_", gm, "-", gu, "_", drug, ".png")
        png(filename, width = 6, height = 4.5, res = 600, units = "in")
        print(p)
        dev.off()

        # perform t-tests
        # todo: save results to dataframe per pset

    }

}

###########################################################
# ttest TODO
###########################################################



dko <- sen$Lapatinib[rownames(sen) %in% rownames(sub[sub$status == "DKO",])]
sko <- sen$Lapatinib[rownames(sen) %in% rownames(sub[sub$status != "DKO",])]
t.test(dko, sko, alternative = "two.sided")

