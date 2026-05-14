#' Plot number of patients
#' 
plot_num_patients <- function(AA_results) {

    toPlot <- table(unique(AA_results[,1:2])$cohort) |> as.data.frame()
    toPlot$Centre <- ifelse(toPlot$Var1 %in% PM2C_cohorts, "PM2C", "BCCA")
    toPlot$Centre <- factor(toPlot$Centre, levels = names(centre_pal))
    toPlot$Var1 <- factor(toPlot$Var1, levels = c(PM2C_cohorts, BCCA_cohorts))

    p <- ggplot(toPlot, aes(x = Var1, y = Freq, fill = Centre)) +
        geom_bar(stat = "identity", color = "black") +
        geom_text(
            aes(label = Freq),
            vjust = -0.3,
            color = "black",
            size = 3
        ) +
        scale_fill_manual(values = centre_pal) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(y = "Number of Smples (with Amplicons)", x = "Cohort")

    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/num_patients.png")
    png(filename, width = 13, height = 4, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Plot number of amplicons in each sample
#' 
plot_num_amplicons <- function(toPlot) {

    toPlot$Centre <- factor(toPlot$Centre, levels = names(centre_pal))
    toPlot$cohort <- factor(toPlot$cohort, levels = c(PM2C_cohorts, BCCA_cohorts))

    p <- ggplot(toPlot, aes(x = cohort, y = Freq, fill = Centre)) +
        geom_boxplot() + 
        #geom_jitter(width = 0.25, size = 0.8, alpha = 0.5) +
        scale_fill_manual(values = centre_pal) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Cohort", y = "Number of Amplicons in Sample")

    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/num_amplicons.png")
    ggsave(filename, p, width = 13, height = 4)

    p <- ggplot(toPlot, aes(x = cohort, y = Freq, fill = Centre)) +
        geom_boxplot(outlier.shape = NA) + 
        scale_y_continuous(limits = c(0, 65)) +
        #geom_jitter(width = 0.25, size = 0.8, alpha = 0.5) +
        scale_fill_manual(values = centre_pal) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Cohort", y = "Number of Amplicons in Sample")

    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/num_amplicons_rmoutliers.png")
    ggsave(filename, p, width = 13, height = 4)
}

plot_num_amplicons2 <- function(num_amplicons) {

    lim <- max(num_amplicons$Freq) + 1

    p <- ggplot(num_amplicons, aes(x = Freq, y = cohort)) +
        geom_density_ridges(
            fill = binary_pal2[1],
            alpha = 0.9,
            scale = 0.9
        ) +
        xlim(0, lim) +
        theme_minimal() +
        labs(x = "Number of Amplicons in Sample", y = "Cohort")

    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/num_amplicons_density.png")
    png(filename, width = 5, height = 6, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Plot number of amplified genes in each sample
#' 
#' @param gene string. "oncogene" or "all"
#' 
plot_num_genes <- function(num_amplicons, gene) {

    col <- switch(
        gene,
        oncogene = "n_genes",
        all = "n_genes_all"
    )
    lab <- switch(
        gene,
        oncogene = "Oncogenes",
        all = "Gene"
    )

    p <- ggplot(num_amplicons, aes(x = .data[[col]], y = cohort)) +
        geom_density_ridges(
            fill = binary_pal2[1],
            alpha = 0.9,
            scale = 0.9
        ) +
        theme_minimal() +
        labs(x = paste("Number of", lab, "in Sample"), y = "Cohort")

    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/num_genes_", lab, ".png")
    png(filename, width = 5, height = 6, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Plot proportion of ecDNA counts across cohorts
#' 
#' @param type string. ecDNA classification ("all" or "binary")
#' 
plot_ecDNA_counts <- function(ecDNA_counts, type = "all") {

    # binarize ecDNA presence
    if (type == "binary") {
        ecDNA_counts$count[ecDNA_counts$count > 0] <- 1
    }

    # plot distribution of all counts
    toPlot <- ecDNA_counts %>%
        count(cohort, count) %>%
        group_by(cohort) %>%
        mutate(prop = n / sum(n) * 100)
    toPlot$count <- factor(toPlot$count, levels = c(max(toPlot$count):0))
    
    # add labels
    if (type == "all") {
        toPlot$label <- ifelse(
            toPlot$count == 0,
            paste0(toPlot$n, "\n(", round(toPlot$prop, 2), "%)"),
            ""
        )
    } else {
        toPlot$label <- paste0(toPlot$n, "\n(", round(toPlot$prop, 2), "%)")
    }
    fill_lab <- ifelse(type == "all", "ecDNA count\nin sample", "ecDNA\nDetection")
    

    # create palette
    pal <- switch(type, all = binary_pal, binary = binary_pal2)
    n <- max(as.integer(as.character(toPlot$count))) + 1
    pal <- colorRampPalette(pal)(n)

    p <- ggplot(toPlot, aes(x = cohort, y = prop, fill = factor(count))) +
        geom_col(color = "black") +
        geom_text(
            aes(label = label),
            position = position_stack(vjust = 0.5),
            color = "black",
            size = 3
        ) +
        scale_fill_manual(values = pal) + 
        theme_classic() +
        labs(x = "Cohort", y = "Proportion of samples (%)", fill = fill_lab)


    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/ecDNA_counts_", type, ".png")
    png(filename, width = 11, height = 4, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Plot proportion of amplicons with oncogenes across cohorts
#' 
plot_oncogene_counts <- function(toPlot) {

    # plot distribution of all counts
    toPlot <- toPlot %>%
        count(cohort, n_genes) %>%
        group_by(cohort) %>%
        mutate(prop = n / sum(n) * 100)
    toPlot$n_genes <- factor(toPlot$n_genes, levels = c(max(toPlot$n_genes):0))
    
    # add labels
    toPlot$label <- ifelse(
        toPlot$n_genes == 0,
        paste0(toPlot$n, "\n(", round(toPlot$prop, 2), "%)"),
        ""
    )

    # create palette
    n <- max(as.integer(as.character(toPlot$n_genes))) + 1
    pal <- colorRampPalette(binary_pal)(n)

    p <- ggplot(toPlot, aes(x = cohort, y = prop, fill = factor(n_genes))) +
        geom_col(color = "black") +
        geom_text(
            aes(label = label),
            position = position_stack(vjust = 0.5),
            color = "black",
            size = 3
        ) +
        scale_fill_manual(values = pal) + 
        theme_classic() +
        labs(x = "Cohort", y = "Proportion of samples (%)", fill = "Number of\nAmplified\nGenes")


    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/oncogene_counts.png")
    png(filename, width = 11, height = 4, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Plot proportion of amplicon classes across cohorts
#' 
plot_amplicon_class <- function(toPlot) {
    toPlot <- toPlot %>%
        group_by(cohort, Classification) %>%
        summarise(total_genes = sum(n_genes), .groups = "drop") %>%
        group_by(cohort) %>%
        mutate(prop = total_genes / sum(total_genes) * 100)
    toPlot$Classification <- factor(toPlot$Classification, levels = names(amplicon_class_pal))

    p <- ggplot(toPlot, aes(x = cohort, y = prop, fill = Classification)) +
            geom_col(color = "black") +
            geom_text(
                aes(label = total_genes),
                position = position_stack(vjust = 0.5),
                color = "black",
                size = 3
            ) +
            scale_fill_manual(values = amplicon_class_pal) + 
            theme_classic() +
            labs(x = "Cohort", y = "Proportion of samples (%)", fill = "Amplicon\nClass")

    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/amplicon_class_amplified_oncogenes.png")
    png(filename, width = 11, height = 4, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Plot proportion of ecDNA amplicons counts across cohorts
#' 
#' @param gene string. "oncogene" or "all"
#' 
plot_ecDNA_genes <- function(ecDNA_amplicons, gene) {

    col <- switch(
        gene,
        oncogene = "Oncogenes",
        all = "All genes"
    )

    ecDNA_amplicons$amplification <- ifelse(
        ecDNA_amplicons[[col]] == "[]",
        0, 1
    )

    # plot distribution of all counts
    toPlot <- ecDNA_amplicons %>%
        count(cohort, amplification) %>%
        group_by(cohort) %>%
        mutate(prop = n / sum(n) * 100)
    toPlot$amplification <- factor(toPlot$amplification, levels = c(1,0))
    toPlot$label <- paste0(toPlot$n, "\n(", round(toPlot$prop, 2), "%)")

    p <- ggplot(toPlot, aes(x = cohort, y = prop, fill = factor(amplification))) +
        geom_col(color = "black") +
        geom_text(
            aes(label = label),
            position = position_stack(vjust = 0.5),
            color = "black",
            size = 3
        ) +
        scale_fill_manual(values = binary_pal2) + 
        theme_classic() +
        labs(x = "Cohort", y = "Proportion of amplicons (%)", fill = "Presence of\nOncogene on\necDNA Amplicon")


    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/genes_", gene, "_on_ecDNA.png")
    png(filename, width = 11, height = 4, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Plot heatmap of oncogenes
#' 
plot_oncogene_heatmap <- function(toPlot, label = "ecDNA") {

    # define colour palettes
    pal <- brewer.pal(brewer.pal.info["Paired", "maxcolors"], "Paired")
    n_cohorts <- length(unique(sub("_.*", "", colnames(toPlot))))
    pal <- pal[1:n_cohorts]
    names(pal) <- unique(sub("_.*", "", colnames(toPlot)))

    # plot heatmap of oncogene detection
    row_ha <- rowAnnotation(num_features = rowSums(toPlot))
    col_ha <- HeatmapAnnotation(
        sample = sub("_amplicon.*", "", colnames(toPlot)),
        cohort = sub("_.*", "", colnames(toPlot)),
        col = list(cohort = pal)
    )

    ht <- Heatmap(
        toPlot, name = "Oncogene\nDetection", 
        top_annotation = col_ha, right_annotation = row_ha,
        show_column_names = FALSE,
        col = c("0" = binary_pal[2], "1" = binary_pal[1]),
        cluster_columns = FALSE,
        rect_gp = gpar(col = "grey80", lwd = 0.5)
    )

    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/", label, "_oncogenes_heatmap.png")
    cat("Saving figure to", filename, "\n")
    png(filename, width = 13, height = 15, res = 600, units = "in")
    print(draw(ht))
    dev.off()
}

#' Plot class comparison of heatmap of oncogenes
#' 
plot_class_oncogene_heatmap <- function(toPlot) {

    ht <- Heatmap(
        toPlot, name = "Oncogene\nDetection",
        rect_gp = gpar(col = "grey80", lwd = 0.5),
        col = magma(256),
        column_names_gp = gpar(fontsize = 7)
    )

    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/class_oncogenes_heatmap.png")
    cat("Saving figure to", filename, "\n")
    png(filename, width = 15, height = 4, res = 600, units = "in")
    print(draw(ht))
    dev.off()
}