# load libraries
suppressPackageStartupMessages({
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
# Validate SKOs in each PSet
###########################################################

validate_SKOs <- function(rna, sen, stat, pset) {

    # load in rna
    filepath <- paste0("data/procdata/psets/", rna, "_RSEM_TPM_reheadered.tsv")
    rna <- fread(filepath, data.table = FALSE)
    rownames(rna) <- rna$Hugo_Symbol
    rna$Hugo_Symbol <- NULL

    # initiate dataframe
    results <- data.frame(matrix(nrow=0, ncol=9))
    colnames(results) <- c(
        "PSet", "Gene_m", "Gene_u", "Drug", "Stat",
        "p.sko1", "p.sko2", "p.dko1", "p.dko2"
    )

    # get stat_dir
    stat_dir <- stat
    if (stat == "1st Qu.") stat_dir <- "quartile1"

    for (i in 1:nrow(top_SKO_tar)) {

        print(i)

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

        # get thirds (non-parametric)
        thirds <- apply(sub, 2, function(x) quantile(x, probs = 1/3, na.rm = TRUE))
        stats$thirds <- c(thirds[[gm]], thirds[[gu]])

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
        sub$status <- factor(sub$status, levels = names(SL_pal))

        # make dataframe of common cells
        common_cells <- intersect(rownames(sub), rownames(sen_targ))
        df <- data.frame(
            sample = common_cells,
            SL = sub$status[match(common_cells, rownames(sub))]
        )
        if (length(df$SL[df$SL == "SKO_targetable"]) < 2) next

        # iterate for every available drug
        for (drug in colnames(sen_targ)) {
            print(paste("--------", drug))
            df$sensitivity <- sen_targ[[drug]][match(common_cells, rownames(sen_targ))]
            toPlot <- df[complete.cases(df),]
            if (length(toPlot$SL[toPlot$SL == "SKO_targetable"]) < 2) next

            # plot boxplots
            p <- ggplot(toPlot, aes(x = SL, y = sensitivity, fill = SL)) +
            geom_boxplot() + geom_jitter(width = 0.2, alpha = 0.5) +
            scale_fill_manual("Synthetic Lethal\nInteraction", values = SL_pal) +
            scale_x_discrete(labels = c("SKO\n(Targeted)", "SKO\n(Other)", "DKO", "no KO")) +
            theme_minimal() +
            theme(panel.border = element_rect()) +
            labs(y = paste(drug, "Response (AAC)"), x = "Synthetic Lethal Interaction") +
            ggtitle(paste0("Dataset: ", pset, ", SKO: ", gm, " & ", gu))
            
            filename <- paste0("data/results/figures/CCL-top-SKOs/", stat_dir, "/", pset, "/", gm, "-", gu, "_", drug, ".png")
            png(filename, width = 5, height = 4, res = 600, units = "in")
            print(p)
            dev.off()

            # get SL classes
            sko_t <- toPlot$sensitivity[toPlot$SL == "SKO_targetable"]
            sko_o <- toPlot$sensitivity[toPlot$SL == "SKO_other"]
            dko <- toPlot$sensitivity[toPlot$SL == "DKO"]
            noko <- toPlot$sensitivity[toPlot$SL == "noKO"]

            # perform t-tests
            sko1 <- ifelse(
                length(sko_t) > 3 & length(sko_o) > 3,
                t.test(sko_t, sko_o, alternative = "greater")$p.value,
                NA
            )
            sko2 <- ifelse(
                length(sko_t) > 3 & length(noko) > 3,
                t.test(sko_t, noko, alternative = "greater")$p.value,
                NA
            )
            dko1 <- ifelse(
                length(dko) > 3 & length(noko) > 3,
                t.test(dko, noko, alternative = "greater")$p.value,
                NA
            )
            dko2 <- ifelse(
                length(dko) > 3 & length(sko_o) > 3,
                t.test(dko, sko_o, alternative = "greater")$p.value,
                NA
            )
        
            row <- data.frame(
                PSet = pset, Gene_m = gm, Gene_u = gu, 
                Drug = drug, Stat = stat, 
                p.sko1 = sko1, p.sko2 = sko2, 
                p.dko1 = dko1, p.dko2 = dko2
            )
            results <- rbind(results, row)
        }
    }
    results <- results[order(results$p.sko1, results$p.dko1),]
    filepath <- paste0("data/results/data/CCL-top-SKOs/", stat_dir, "/", pset, ".csv")
    write.csv(results, file = filepath, quote = FALSE, row.names = FALSE)

    return(results)
}

stat_options <- c("1st Qu.", "Median", "thirds")

s = 3
ccle <- validate_SKOs("CCLE", ctrp_sen, stat_options[s], "CTRP")
gdsc <- validate_SKOs("GDSC", gdsc_sen, stat_options[s], "GDSC")
gcsi <- validate_SKOs("gCSI", gcsi_sen, stat_options[s], "gCSI")
gray <- validate_SKOs("GRAY", gray_sen, stat_options[s], "GRAY")
