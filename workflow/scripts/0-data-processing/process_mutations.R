# load libraries
suppressPackageStartupMessages({
    library(maftools)
    library(EnsDb.Hsapiens.v86)
    library(AnnotationDbi)
})

setwd("/cluster/projects/bhklab/projects/SL_MOHCCN/data/rawdata/ecDNA/mafs")

###########################################################
# Make mutation counts matrix
###########################################################

# get mafs
message("loading mafs")
mafs <- merge_mafs(list.files())

# create counts matrix
message("creating matrix")
mut <- mutCountMatrix(
  mafs,
  includeSyn = FALSE,
  countOnly = NULL,
  removeNonMutated = FALSE
  #removeNonMutated = TRUE
)

# save mutation counts matrix
message("saving matrix")
write.table(
    mut, 
    file = "/cluster/projects/bhklab/projects/SL_MOHCCN/data/procdata/mutations/mat.tsv", 
    quote = F, sep = "\t", col.names = T, row.names = T
)

###########################################################
# Get consensus
###########################################################

setwd("/cluster/projects/bhklab/projects/SL_MOHCCN/data/procdata/mutations")

files <- list.files(pattern = "tsv")
consensus <- c()
for (file in files) {
    x <- read.table(file)
    consensus <- unique(c(consensus, rownames(x)))
}
saveRDS(consensus, file = "final/consensus.RDS")


###########################################################
# Make final matrix
###########################################################

consensus <- readRDS("final/consensus.RDS")

mut_df <- data.frame(matrix(nrow=0, ncol=length(consensus)))

# merge matrices
for (file in files) {
    mut <- read.table(file)
    mut <- as.data.frame(t(mut))
    missing_genes <- consensus[-which(consensus %in% colnames(mut))]
    mut[missing_genes] <- 0
    mut <- mut[,match(consensus, colnames(mut))]
    mut_df <- rbind(mut_df, mut)
}

write.table(
    mut_df, 
    file = "/cluster/projects/bhklab/projects/SL_MOHCCN/data/procdata/mutations/final/PM_mutations.tsv", 
    quote = F, sep = "\t", col.names = T, row.names = T
)

###########################################################
# Make TP53 mutation counts matrix
###########################################################

setwd("/cluster/projects/bhklab/projects/SL_MOHCCN/data/procdata/mutations")

df <- data.frame(matrix(nrow=1, ncol=0))
for (file in list.files()) {
    x <- read.table(file)
    if ("TP53" %in% rownames(x)) {
        x <- x["TP53",]
        df <- cbind(df, x)
    }
}
df <- as.data.frame(t(df))
df$sample <- rownames(df)
write.table(
    df,
    file = "/cluster/projects/bhklab/projects/SL_MOHCCN/data/procdata/mutations/TP53.tsv", 
    quote = F, sep = "\t", col.names = T, row.names = T
)

###########################################################
# Misc
###########################################################

#mat1: 1:50
#mat2: 51:100
#mat3: 101:150
#mat4: 151:200
#mat5: 201:250
#mat6: 251:300
#mat7: 301:350
#mat8: 351:400
#mat9: 401:450
#mat10: 451:500
#mat11: 501:550
#mat12: 551:600
#mat13: 601:615
#mat38: 616:630
#mat35: 631:650
#mat14: 651:670
#mat34: 671:700
#mat15: 701:750
#mat16: 751:800
#mat17: 801:850
#mat18: 851:870
#mat37: 871:900
#mat19: 901:950
#mat20: 951:1000
#mat21: 1001:1015
#mat42: 1016:1025
#mat39: 1026:1035
#mat44: 1036:1050
#mat22: 1051:1070
#mat40: 1071:1077
#mat46: 1078:1085
#mat41: 1086:1100
#mat23: 1101:1115
#mat47: 1116:1125
#mat43: 1126:1135
#mat48: 1136:1143
#mat51: 1144:1150
#mat24: 1151:1162
#mat52: 1163:1175
#mat45: 1176:1180
#mat54: 1181:1185
#mat53: 1186:1200
#mat25: 1201:1250
#mat26: 1251:1300
#mat27: 1301:1350
#mat28: 1351:1400
#mat29: 1401:1450
#mat30: 1451:1475
#mat49: 1476:1485
#mat55: 1486:1500
#mat31: 1501:1514
#mat56: 1515:1525
#mat50: 1526:1550
#mat32: 1551:1555
#mat57: 1556:1575
#mat58: 1576:1585
#mat59: 1586:1600
#mat33: 1601:1650
#mat36: 1651:1689


cols <- readRDS("data/procdata/mutations/consensus.RDS")
enst_ids <- grep("^ENST", cols, value = TRUE)
gene_symbols <- setdiff(cols, enst_ids)

# strip version suffixes if present, e.g. ENST00000424460.5 -> ENST00000424460
enst_clean <- sub("\\..*", "", enst_ids)

mapping <- select(EnsDb.Hsapiens.v86, keys = enst_clean, keytype = "TXID", columns = "SYMBOL")
colnames(mapping) <- c("ensembl_transcript_id", "hgnc_symbol")

overlap <- mapping[mapping$hgnc_symbol %in% gene_symbols, ]

enst_to_symbol <- overlap$hgnc_symbol
names(enst_to_symbol) <- overlap$ensembl_transcript_id
saveRDS(enst_to_symbol, file = "data/procdata/mutations/enst_to_symbol.RDS")
