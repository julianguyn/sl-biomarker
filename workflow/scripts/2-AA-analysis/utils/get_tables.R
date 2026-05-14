#' Function to process cohort names from input files
#' 
#' Requirements: 
#' source("workflow/scripts/utils/cohorts.R")
#' Column 'cohort' from parent directory of AA output folders
#' 
process_cohort_names <- function(df) {

  PM2C_results <- df[df$cohort %in% PM2C_cohorts,]
  BCCA_results <- df[-which(df$cohort %in% PM2C_cohorts),]

  for (cohort in BCCA_cohorts) {
    BCCA_results$cohort[grep(cohort, BCCA_results$cohort)] <- cohort
  }

  PM2C_results$Centre <- "PM2C"
  BCCA_results$Centre <- "BCCA"

  df <- rbind(PM2C_results, BCCA_results)

  # factor for plotting
  df$Centre <- factor(df$Centre, levels = c("PM2C", "BCCA"))
  df$cohort <- factor(df$cohort, levels = c(PM2C_cohorts, BCCA_cohorts))

  return(df)
}