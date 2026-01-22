# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(stringr)
    library(purrr)
    library(dplyr)
})

RAWDATA_DIR <- "/cluster/projects/bhklab/projects/SL_MOHCCN/data/rawdata/ecDNA/PM2C_Batch_1/"
PROCDATA_DIR <- "/cluster/projects/bhklab/projects/SL_MOHCCN/data/procdata/"

###########################################################
# Load in data
###########################################################

# get samples (directory names)
sample_names <- list.files(paste0(RAWDATA_DIR, "Amplicon_Architect"))
samples <- data.frame(
    cohort = sub("_.*", "", sample_names),
    sample = sample_names
)

## TODO: missing BTC files
samples <- samples[samples$cohort != "BTC", ]

###########################################################
# Get AA summary text files
###########################################################

# get AA summary text files
samples$summary <- paste0(
    RAWDATA_DIR, "Amplicon_Architect/",
    samples$sample, "/",
    samples$sample, "_AA_results/",
    samples$sample, "_summary.txt"
)

# check all files exist
#all_files <- list.files(
#    paste0(RAWDATA_DIR, "Amplicon_Architect"),
#    recursive = TRUE, 
#    full.names = TRUE,
#    pattern = "*_summary.txt"
#)
#checks <- samples$summary %in% all_files
#if (all(checks) != TRUE) {
#    message("The following summary files do not exist:")
#    print(samples$summary[which(checks == FALSE)])
#}

###########################################################
# Get classification data filepaths
###########################################################

# files to get
files <- c(
    #"amp_classification" = "amplicon_classification_profiles.tsv",
    "result_table" = "result_table.tsv",
    #"gene_list" = "gene_list.tsv",
    #"feature_basic_properties" = "feature_basic_properties.tsv",
    #"feature_entropy" = "feature_entropy.tsv",
    "ecDNA_counts" = "ecDNA_counts.tsv"#,
    #"ecDNA_context_calls" = "ecDNA_context_calls.tsv"
)

# build dataframe of filepaths
missing <- c()
for (file in names(files)) {

    # make column
    samples[[file]] <- NA_character_

    for (i in 1:nrow(samples)) {
        # get filepath
        samples[[file]][i] <- paste0(
            RAWDATA_DIR, "Amplicon_Architect/",
            samples$sample[i], "/",
            samples$sample[i], "_classification/",
            samples$sample[i], "_", files[file]
        )
    }

    # check all files exist
    all_files <- list.files(
        paste0(RAWDATA_DIR, "Amplicon_Architect"),
        recursive = TRUE,
        full.names = TRUE,
        pattern = paste0("*", files[file])
    )
    checks <- samples[[file]] %in% all_files
    if (all(checks) != TRUE) {
        message("The following output files do not exist:")
        f <- samples[[file]][which(checks == FALSE)]
        print(sub("/cluster/projects/bhklab/projects/SL_MOHCCN/data/rawdata/ecDNA/PM2C_Batch_1/Amplicon_Architect/", "", f))
        missing <- c(missing, f)
    }
}

###########################################################
# Compile AA results
###########################################################

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
ecDNA_amplicons <- AA_results[Classification == "ecDNA"]

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