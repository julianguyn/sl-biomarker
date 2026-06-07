# process rsem on H4H
library(data.table)

dirs <- c(
    "../../data/rawdata/ecDNA/PM2C_Batch_1/RSEM/",
    "../../data/rawdata/ecDNA/PM2C_Batch_2/RSEM/"
)
counts <- data.frame(matrix(nrow=60603, ncol=0))
get_genes <- FALSE

for (dir in dirs) {
    files <- list.files(dir, full.names = TRUE)
    print(length(files))

    for (file in files) {
        sample <- sub(paste0(dir, "/"), "", sub(".RSEM.genes.results", "", file))
        df <- fread(file)
        counts[[sample]] <- df$TPM
        if (get_genes == FALSE) {
            rownames(counts) <- df$gene_id
            get_genes <- TRUE
        }
    }
}

saveRDS(counts, file = "rsem.rds") #1534 samples

# save to data/procdata/RNA_TPM/rsem.rds
# process rsem


PM_meta <- read.csv("metadata/2026-05-29_PM_meta.csv")
#BC_meta <- read.csv("metadata/2026-05-29_BC_meta.csv")

pm_rna <- readRDS("data/procdata/RNA_TPM/rsem.rds")
pm_rna <- pm_rna[,colnames(pm_rna) %in% PM_meta$RNA] #1521 samples
colnames(pm_rna) <- PM_meta$Sample_ID[match(colnames(pm_rna), PM_meta$RNA)]

saveRDS(pm_rna, file = "data/procdata/RNA_TPM/PM_RNA.rds")
