# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(purrr)
    library(dplyr)
})

###########################################################
# Set variables
###########################################################

# data directories
RAWDATA_DIR <- "/cluster/projects/bhklab/projects/SL_MOHCCN/data/rawdata/ecDNA/"
PROCDATA_DIR <- "/cluster/projects/bhklab/projects/SL_MOHCCN/data/procdata/"

# batches
batches <- c("PM2C_Batch_1", "PM2C_Batch_2", "BCCA_Batch_1", "BCCA_Batch_2")
# 461, 1275, 781, 1512
# total: 4029 samples

###########################################################
# Load in data
###########################################################

# get samples (directory names)
all_sample_names <- c()
samples <- data.frame(matrix(nrow=0, ncol=4))
for (batch in batches) {
    if (batch != "BCCA_Batch_2") {
        sample_names <- list.files(paste0(RAWDATA_DIR, batch, "/Amplicon_Architect"))
        RESULT_DIR <- "/Amplicon_Architect/"
    } else {
        sample_names <- list.files(paste0(RAWDATA_DIR, batch, "/AA_OUTPUT"))
        RESULT_DIR <- "/AA_OUTPUT/"
    }
    all_sample_names <- c(all_sample_names, sample_names)
    samples <- rbind(samples, data.frame(
        batch = batch,
        rdir = RESULT_DIR,
        cohort = sub("_.*", "", sample_names),
        sample = sample_names
    ))
}

## TODO: missing BTC files
#samples <- samples[samples$cohort != "BTC", ]

###########################################################
# Manually extract files
###########################################################

samples$dir <- paste0(RAWDATA_DIR, samples$batch, samples$rdir)
samples$result_table <- NA_character_
samples$ecDNA_counts <- NA_character_

for (i in 1:nrow(samples)) {

    # result_table
    result_table <- list.files(
        paste0(samples$dir[i], samples$sample[i], "/", samples$sample[i], "_classification/"),
        recursive = TRUE, full.names = TRUE,
        pattern = paste0("result_table.tsv"))
    samples$result_table[i] <- ifelse(length(result_table) > 0, result_table, "missing")

    # ecDNA_counts
    ecDNA_counts <- list.files(
        paste0(samples$dir[i], samples$sample[i], "/", samples$sample[i], "_classification/"),
        recursive = TRUE, full.names = TRUE,
        pattern = paste0("ecDNA_counts.tsv"))
    samples$ecDNA_counts[i] <- ifelse(length(ecDNA_counts) > 0, ecDNA_counts, "missing")
}
saveRDS(samples, file = paste0(PROCDATA_DIR, "AA_output/sample_df.rds"))

###########################################################
# Compile AA results
###########################################################

samples <- samples[samples$result_table != "missing",]

# rbind all result table files
message("\nCompiling AA results")
AA_results <- samples$result_table %>%
  map_dfr(fread) %>%
  mutate(
    cohort = sub("_.*", "", 'Sample name'),
    .before = 1
  )

message("Saving AA_results")
saveRDS(
    AA_results,
    file = paste0(PROCDATA_DIR, "AA_output/AA_results.rds")
)

# get ecDNA-classified amplicons
ecDNA_amplicons <- AA_results[AA_results$Classification == "ecDNA"]

message("Saving ecDNA AA_results")
saveRDS(
    ecDNA_amplicons,
    file = paste0(PROCDATA_DIR, "AA_output/ecDNA_amplicons.rds")
)

###########################################################
# Get ecDNA counts
###########################################################

# initiate dataframe
ecDNA_counts <- data.frame(
    cohort = samples$cohort,
    sample = samples$sample,
    count = NA
)

# get ecDNA counts
for (i in 1:length(samples$ecDNA_counts)) {
    counts <- fread(samples$ecDNA_counts[i])
    ecDNA_counts$count[i] <- counts$ecDNA_count
}

message("Saving ecDNA_counts")
write.csv(
    ecDNA_counts,
    file = paste0(PROCDATA_DIR, "AA_output/ecDNA_counts.csv"),
    quote = FALSE,
    row.names = FALSE
)