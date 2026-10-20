# table S2 and figure S2

###########################################################
# Process metadata
###########################################################

# load PM oncotree
mohda <- read.csv("metadata/mohda2026003_oncotree.csv")
mohda$matchID <- sub("WG_", "", mohda$samples)

mohda <- mohda[-grep("Removed", mohda$notes),]

cohort <- table(mohda$cohort) |> as.data.frame()
colnames(cohort) <- c("Cohort", "Number of Patients")

tissue <- table(mohda$oncotree_root) |> as.data.frame()
colnames(tissue) <- c("Oncotree Root", "Number of Patients")

write.csv(cohort, file = "data/results/data/ecDNA_collab/tables2_cohorts.csv", quote = F, row.names = F)
write.csv(tissue, file = "data/results/data/ecDNA_collab/tables2_tissues.csv", quote = F, row.names = F)
