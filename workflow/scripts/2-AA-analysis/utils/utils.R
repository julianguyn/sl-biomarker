#' Helper function to get amplified oncogenes per class
#' 
#' Used in scripts/AA_data_exploration.R
#' @param df data.frame. AA_results from PROCDATA/AA_output/AA_results.rds
#' @param classification string. Type of classification to subset for
#' options include: "BFB", "Complex-non-cyclic", "ecDNA", "Linear"
#' 
get_amplified_oncogenes <- function(df, classification) {

  df <- df[df$Classification == classification,]

  oncogenes <- unlist(
    strsplit(gsub("\\[|\\]|'", "", df$Oncogenes), ",\\s*")
  ) |> unique()

  # initate dataframe to store results
  toPlot <- data.frame(matrix(nrow=nrow(df), ncol=length(oncogenes)))
  colnames(toPlot) <- c(oncogenes)
  rownames(toPlot) <- df$'Feature ID'

  # binarize oncogene detection
  for (gene in oncogenes) {
    toPlot[[gene]][grep(gene, df$Oncogenes)] <- 1
  }

  # format for plotting
  toPlot[is.na(toPlot)] <- 0
  toPlot <- t(toPlot) |> as.data.frame()

  # plot heatmap
  plot_oncogene_heatmap(toPlot, label = classification)

  # get counts of amplified oncogenes
  count_amp <- rowSums(toPlot)
  count_amp <- count_amp[order(count_amp, decreasing = TRUE)]
  count_amp <- data.frame(
    Gene = names(count_amp),
    Count = unname(count_amp)
  )
  colnames(count_amp)[2] <- classification
  return(count_amp)
}
