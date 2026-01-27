suppressPackageStartupMessages({
    library(PharmacoGx)
    library(dplyr)
})

# ------------------------------------------------------
# Helper function to extract PSet gene expression counts

get_pset_SE <- function(pset) {

    pset_options <- c("GRAY", "gCSI", "CCLE", "GDSC2")
    if (!pset %in% pset_options) {
        message("Option not valid. Pick from: GRAY, gCSI, CCLE, GDSC2")
        return(NA)
    }
    if (pset == "GRAY") pset <- "PSet_GRAY2017"
    if (pset == "GDSC2") pset <- "GDSC2-8.2"

    pset_path <- paste0("../BCaATAC/data/rawdata/psets/", pset, ".rds")
    pset <- readRDS(pset_path) |> updateObject()

    # get TPM counts
    rna <- summarizeMolecularProfiles(
        pset, 
        mDataType = "Kallisto_0.46.1.rnaseq"
    ) 
    return(rna)
}

get_gene_meta <- function(rna) {
    keep <- c("seqnames", "start", "end", "strand", "gene_id", "gene_name")
    gene_meta <- rna@elementMetadata[, keep]
    return(gene_meta)
}

get_cell_meta <- function(rna) {
    keep <- c("Cell_line", "sampleid", "batchid", "Tissue_supergroup", "Primary_Tissue", "tissueid")
    keep <- keep[keep %in% colnames(rna@colData)]
    cell_meta <- rna@colData[, keep]
    cell_meta$names <- rownames(rna@colData)
    return(cell_meta)
}

get_pset_rna <- function(rna, gene_meta) {

    counts <- rna |> assay() |> as.data.frame()
    colnames(counts) <- rownames(rna@colData)

    # average across gene symbols
    counts$Genes <- gene_meta$gene_name
    df <- counts %>%
        dplyr::group_by(Genes) %>%
        dplyr::summarise(across(everything(), mean)) %>%
        as.data.frame()
    rownames(df) <- df$Genes
    df$Genes <- NULL

    # undo log transformation
    df <- 2^df - 0.001

    return(df)
}

write_rna <- function(pset_rna, label) {
    df <- cbind(rownames(pset_rna), pset_rna)
    colnames(df)[1] <- "Hugo_Symbol"
    filename <- paste0("data/procdata/psets/", label, "_RSEM_TPM_reheadered.tsv")
    write.table(df, file = filename, quote = FALSE, sep = "\t", row.names = FALSE)
}

###########################################################
# Get TPM count matrices
###########################################################

# get SE objects
gray <- get_pset_SE("GRAY")
gcsi <- get_pset_SE("gCSI")
ccle <- get_pset_SE("CCLE")
gdsc <- get_pset_SE("GDSC2")

# get metadata
gray_gene <- get_gene_meta(gray)
gcsi_gene <- get_gene_meta(gcsi)
ccle_gene <- get_gene_meta(ccle)
gdsc_gene <- get_gene_meta(gdsc)

gray_cell <- get_cell_meta(gray)
gcsi_cell <- get_cell_meta(gcsi)
ccle_cell <- get_cell_meta(ccle)
gdsc_cell <- get_cell_meta(gdsc)

# get TPMs
gray_rna <- get_pset_rna(gray, gray_gene)
gcsi_rna <- get_pset_rna(gcsi, gcsi_gene)
ccle_rna <- get_pset_rna(ccle, ccle_gene)
gdsc_rna <- get_pset_rna(gdsc, gdsc_gene)

# write out gene matrices for SL pipeline
write_rna(gray_rna, "GRAY")
write_rna(gcsi_rna, "gCSI")
write_rna(ccle_rna, "CCLE")
write_rna(gdsc_rna, "GDSC")

# save metadata objects
save(
  gray_gene, gcsi_gene, ccle_gene, gdsc_gene,
  gray_cell, gcsi_cell, ccle_cell, gdsc_cell,
  file = "data/procdata/psets/rna_metadata.RData"
)
