# quick script to look for variables of interest in CANDIG pull

library(data.table)


PROCDATA_DIR <- "data/procdata/"

source("workflow/scripts/2-AA-analysis/utils/utils.R")
source("workflow/scripts/utils/cohorts.R")
source("workflow/scripts/utils/ID_MAPPING.R")

###########################################################
# Load in data
###########################################################

# load in all samples
all_samples <- readRDS(paste0(PROCDATA_DIR, "AA_output/sample_df.rds"))

# process cohort names (from utils/get_tables.R)
all_samples <- process_cohort_names(all_samples)

all_samples$Patient_ID <- sub("^([A-Za-z0-9]+_[0-9]+)_.*$", "\\1", all_samples$sample)
all_samples <- all_samples[,c("Patient_ID", "sample", "cohort", "Centre")]

###########################################################
# Process variables
###########################################################

variables <- data.frame(matrix(nrow=0, ncol=2))
colnames(variables) <- c("file", "variable")


dir <- "data/rawdata/clinical_data_20260923/"
for (file in list.files(dir)) {
    df <- read.csv(paste0(dir, file))
    tobind <- data.frame(file = rep(file, ncol(df)), variable = colnames(df))
    variables <- rbind(variables, tobind)
}


needed_variables <- c(
    "clinical_stage_group",
    "pathological_stage_group",
    "clinical_t_category",
    "clinical_n_category",
    "clinical_m_category",
    "pathological_t_category",
    "pathological_n_category",
    "pathological_m_category",
    "date_of_birth",
    "sex_at_birth",
    "gender",
    "cancer_type_code",
    "primary_site",
    "specimen_anatomic_location",
    "date_of_relapse",
    "date_of_death",
    "date_of_followup",
    "tobacco_smoking_status",
    "tobacco_type",
    "pack_years_smoked"
)

table(needed_variables %in% variables$variable)

variables <- variables[variables$variable %in% needed_variables,]
variables$num_missing <- NA
variables$incomplete_code <- NA

# manual check data completeness
missing_IDs <- c()

get_var <- function(var){
    file <- variables$file[variables$variable == var]

    df <- read.csv(paste0(dir, file))
    var <- df[,colnames(df) %in% c(var, "submitter_donor_id")]
    return(var)
}

na_vars <- c(
    "clinical_stage_group", "pathological_stage_group",
    "sex_at_birth", "gender", 
    "cancer_type_code", "primary_site", "specimen_anatomic_location",
    "tobacco_smoking_status"
)

blank_vars <- c(
    "clinical_t_category", "clinical_n_category",
    "clinical_m_category", "pathological_t_category",
    "pathological_n_category", "pathological_m_category",
    "date_of_birth", "date_of_relapse", "date_of_death",
    "date_of_followup", "tobacco_type"
)

for(var in na_vars){
    x <- get_var(var)
    x <- x[grep("Not available", x[[var]]),]
    variables$num_missing[variables$variable == var] <- nrow(x)
    variables$incomplete_code[variables$variable == var] <- "Not available"

    missing_IDs[[var]] <- x$submitter_donor_id
}

for(var in blank_vars){
    x <- get_var(var)
    x <- x[x[[var]] == "",]
    variables$num_missing[variables$variable == var] <- nrow(x)
    variables$incomplete_code[variables$variable == var] <- ""

    missing_IDs[[var]] <- x$submitter_donor_id
}

write.csv(variables, file = "data/procdata/clinical/clinical_candig_missing_vars.csv", quote = F, row.names = F)
saveRDS(missing_IDs, file = "data/procdata/clinical/clinical_candig_missing_IDs.RDS")

# sanity check
for (var in variables$variable) {
    ids <- length(missing_IDs[[var]])
    print(var)
    print(ids == variables$num_missing[variables$variable == var])
}

# get all unique IDs missing
missing <- c()
for (var in variables$variable) {
    ids <- missing_IDs[[var]]
    missing <- unique(c(missing, ids))
}

# NA for "pack_years_smoked", but minimum is 0.05 -> NA = never smoked?