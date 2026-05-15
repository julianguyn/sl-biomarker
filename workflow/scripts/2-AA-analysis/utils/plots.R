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

    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/exploration/num_patients.png")
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

    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/exploration/num_amplicons.png")
    ggsave(filename, p, width = 13, height = 4)

    p <- ggplot(toPlot, aes(x = cohort, y = Freq, fill = Centre)) +
        geom_boxplot(outlier.shape = NA) + 
        scale_y_continuous(limits = c(0, 120)) +
        #geom_jitter(width = 0.25, size = 0.8, alpha = 0.5) +
        scale_fill_manual(values = centre_pal) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Cohort", y = "Number of Amplicons in Sample")

    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/exploration/num_amplicons_rmoutliers.png")
    ggsave(filename, p, width = 13, height = 4)
}

#' Plot proportion of amplicon classes across cohorts
#' 
#' @param label string. One of: c(all_amplicons, amplified_oncogenes, amplified_all_genes)
plot_amplicon_class <- function(toPlot, label) {

    # amplicon class of amplified genes
    if (label == "all_amplicons") {
        toPlot <- toPlot %>%
            group_by(cohort, Classification) %>%
            summarise(total = n(), .groups = "drop") %>%
            group_by(cohort) %>%
            mutate(prop = total / sum(total) * 100)
        dir <- "exploration"
    } else if (label == "amplified_oncogenes") {
        toPlot <- toPlot %>%
            group_by(cohort, Classification) %>%
            summarise(total = sum(n_genes), .groups = "drop") %>%
            group_by(cohort) %>%
            mutate(prop = total / sum(total) * 100)
        dir <- "gene_analysis"
    } else if (label == "amplified_all_genes") {
        toPlot <- toPlot %>%
            group_by(cohort, Classification) %>%
            summarise(total = sum(n_genes_all), .groups = "drop") %>%
            group_by(cohort) %>%
            mutate(prop = total / sum(total) * 100)
        dir <- "gene_exploration"
    }
    

    toPlot$Classification <- factor(toPlot$Classification, levels = names(amplicon_class_pal))

    p <- ggplot(toPlot, aes(x = cohort, y = prop, fill = Classification)) +
            geom_col(color = "black") +
            geom_text(
                aes(label = total),
                position = position_stack(vjust = 0.5),
                color = "black",
                size = 3
            ) +
            scale_fill_manual(values = amplicon_class_pal) + 
            theme_classic() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = "Cohort", y = "Proportion of samples (%)", fill = "Amplicon\nClass")

    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/", dir, "/amplicon_class_", label, ".png")
    print(filename)
    png(filename, width = 13, height = 4, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Plot distribution of feature median / maximum copy number per amplicon class
#' 
#' @param variable string. 'Feature median copy number' or 'Feature maximum copy number'
#' 
plot_fmcn <- function(toPlot, variable) {

  label <- ifelse(
    variable == 'Feature median copy number',
    "fmedcn", "fmaxcn"
  )

  toPlot$Classification <- factor(toPlot$Classification, levels = names(amplicon_class_pal))

  p <- ggplot(toPlot, aes(x = Classification, y = .data[[variable]], fill = Classification)) +
    geom_boxplot() +
    scale_y_continuous(limits = c(0, 275)) +
    facet_wrap(~cohort, ncol = 10) +
    scale_fill_manual(values = amplicon_class_pal) +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  filename <- paste0(RESULTS_DIR, "figures/AA_exploration/exploration/", label, ".png")
  ggsave(filename, p, w = 15, h = 6)

  p <- ggplot(toPlot, aes(x = Classification, y = .data[[variable]], fill = Classification)) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_continuous(limits = c(0, 80)) +
    facet_wrap(~cohort, ncol = 10) +
    scale_fill_manual(values = amplicon_class_pal) +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  filename <- paste0(RESULTS_DIR, "figures/AA_exploration/exploration/", label, "_no_outliers.png")
  ggsave(filename, p, w = 15, h = 6)
}

#' Plot proportion of samples with ecDNA
#' 
#' @param type string. ecDNA classification ("all" or "binary")
#' 
plot_ecDNA_prop <- function(ecDNA_counts, type = "all") {

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
            toPlot$count == 1,
            paste0(toPlot$n, "\n(", round(toPlot$prop), "%)"),
            ""
        )
    } else {
        toPlot$label <- paste0(toPlot$n, "\n(", round(toPlot$prop), "%)")
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
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Cohort", y = "Proportion of samples (%)", fill = fill_lab)


    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/exploration/ecDNA_counts_", type, ".png")
    png(filename, width = 15, height = 4, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Plot proportion of amplicons with oncogenes across cohorts
#' 
plot_oncogene_counts <- function(toPlot) {

    # plot distribution of all counts
    toPlot$have_oncogene <- ifelse(toPlot$n_genes > 0, 1, 0)
    toPlot <- toPlot %>%
        count(cohort, have_oncogene) %>%
        group_by(cohort) %>%
        mutate(prop = n / sum(n) * 100)
    toPlot$have_oncogene <- factor(toPlot$have_oncogene, levels = c(1,0))
    
    # add labels
    toPlot$label <- paste0(toPlot$n, "\n(", round(toPlot$prop), "%)")

    # create palette
    p <- ggplot(toPlot, aes(x = cohort, y = prop, fill = factor(have_oncogene))) +
        geom_col(color = "black") +
        geom_text(
            aes(label = label),
            position = position_stack(vjust = 0.5),
            color = "black",
            size = 3
        ) +
        scale_fill_manual(values = binary_pal2) + 
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Cohort", y = "Proportion of samples (%)", fill = "Amplified\nOncogene(s)")


    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/gene_exploration/prop_have_oncogene.png")
    png(filename, width = 15, height = 4, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Plot distribution of feature median/maximum CN across amplified oncogenes
#' by amplicon class
#' 
#' @param tukey test ouputs for median and maximum CN by amplicon class
#' 
plot_fmcn_genes <- function(toPlot, med_tukey, max_tukey, label) {

    toPlot$Classification <- factor(toPlot$Classification, levels = names(amplicon_class_pal))

    p1 <- ggplot(toPlot, aes(x = Classification, y = log2(Feature_median_copy_number), fill = Classification)) +
        geom_violin() + geom_boxplot(width = 0.15) +
        scale_fill_manual(values = amplicon_class_pal) +
        theme_bw() +
        theme(legend.position = "none") +
        labs(y = paste0("log2(Feature median copy number)"), title = "Feature median copy number")

    p2 <- ggplot(toPlot, aes(x = Classification, y = log2(Feature_maximum_copy_number), fill = Classification)) +
        geom_violin() + geom_boxplot(width = 0.15) +
        scale_fill_manual(values = amplicon_class_pal) +
        theme_bw() +
        theme(legend.position = "none") +
        labs(y = paste0("log2(Feature maximum copy number)"), title = "Feature maximum copy number")

    if (label != "all_genes") {
        p1 <- p1 + geom_signif(
                y_position = c(6.45, 6.1, 8.5, 8, 7.5),
                xmin = c(1, 2, 1, 2, 3), xmax = c(2, 3, 4, 4, 4),
                annotation = med_tukey$'p adj',
                tip_length = 0.01,
                textsize = 4
        )
        p2 <- p2 + geom_signif(
                y_position = c(9.25, 9.75, 8.55, 8.15, 9.05),
                xmin = c(1, 1, 2, 2, 3), xmax = c(2, 3, 4, 3, 4),
                annotation = max_tukey$'p adj',
                tip_length = 0.01,
                textsize = 4
        )
    }

    p <- p1 + p2
    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/gene_exploration/amplified_", label, "_CN_amplicon_class.png")
    ggsave(filename, p, w = 10, h = 6)
}

#' Plot proportion of ecDNA amplicons counts across cohorts
#' 
#' @param gene string. "oncogene" or "all"
#' 
plot_ecDNA_genes <- function(ecDNA_amplicons, gene) {

    col <- switch(
        gene,
        oncogene = "Oncogenes",
        all = "All_genes"
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
    toPlot$label <- paste0(toPlot$n, "\n(", round(toPlot$prop), "%)")

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
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Cohort", y = "Proportion of amplicons (%)",
             fill = paste0("Presence of\n", col, " on\necDNA Amplicon"))


    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/gene_exploration/genes_", gene, "_on_ecDNA.png")
    png(filename, width = 13, height = 4, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Plot top oncogenes across amplicon classes and CN
#' 
plot_top_amp_oncogenes <- function(top_oncogene, oncogene_df, sig_ecdna_linear) {

    onco_order <- top_oncogene$Var1
    top_oncogene$Var1 <- factor(top_oncogene$Var1, levels = onco_order)

    # plot count of oncogenes
    p1 <- ggplot(top_oncogene, aes(x = Var1, y = Freq)) +
        geom_bar(stat = "identity") +
        geom_text(
                aes(label = Freq),
                vjust = -0.3,
                color = "black",
                size = 3
            ) +
        scale_y_continuous(limits = c(0, 73)) +
        theme_bw() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank()
        ) +
        labs(y = "No. Amplicons\nw Oncogene")

    # plot proportion of amplicon class with oncogene
    toPlot <- oncogene_df[oncogene_df$oncogene %in% onco_order,]
    toPlot <- toPlot %>%
        group_by(oncogene, Classification) %>%
        summarise(total = n(), .groups = "drop") %>%
        group_by(oncogene) %>%
        mutate(prop = total / sum(total) * 100)
    toPlot$oncogene <- factor(toPlot$oncogene, levels = onco_order)
    toPlot$Classification <- factor(toPlot$Classification, levels = names(amplicon_class_pal))
    p2 <- ggplot(toPlot, aes(x = oncogene, y = prop, fill = Classification)) +
        geom_col(color = "black") +
        geom_text(
            aes(label = total),
            position = position_stack(vjust = 0.5),
            color = "black",
            size = 3
        ) +
        scale_fill_manual(values = amplicon_class_pal) + 
        theme_bw() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank()
        ) +
        labs(x = "Oncogene", y = "Proportion of\namplicons (%)", fill = "Amplicon\nClass")

    # heatmap of sig_ecdna_linear
    toPlot <- data.frame(oncogene = onco_order, sig = ifelse(onco_order %in% sig_ecdna_linear, 1, 0))
    toPlot$oncogene <- factor(toPlot$oncogene, levels = onco_order)
    toPlot$sig <- factor(toPlot$sig, levels = c(1:0))

    p3 <- ggplot(toPlot, aes(x = oncogene, y = "", fill = sig)) +
        geom_tile(color = "black") +
        scale_fill_manual(values = c(unname(amplicon_class_pal["ecDNA"]), "white")) +
        theme_void() +
        theme(legend.position = "none")

    # helper function to plot distribution of feature median CN
    helper_plot_fmcn_class <- function(df, class) {

        toPlot <- df[df$oncogene %in% onco_order,]
        toPlot <- toPlot[toPlot$Classification == class,]
        toPlot$oncogene <- factor(toPlot$oncogene, levels = onco_order)

        p <- ggplot(toPlot, aes(x = oncogene, y = log2(Feature_median_copy_number))) +
            geom_boxplot(fill = amplicon_class_pal[class]) +
            geom_jitter(width = 0.2, alpha = 0.3) +
            scale_x_discrete(drop = FALSE) +
            scale_y_continuous(limits = c(2, 7.5)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
            labs(y = "log2(Median\nCopy Number)", x = "Oncogene")
        if (class != "ecDNA") {
            p <- p + theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank()
            )
        }
        return(p)
    }

    p4 <- helper_plot_fmcn_class(oncogene_df, "BFB")
    p5 <- helper_plot_fmcn_class(oncogene_df, "Complex-non-cyclic")
    p6 <- helper_plot_fmcn_class(oncogene_df, "Linear")
    p7 <- helper_plot_fmcn_class(oncogene_df, "ecDNA")

    p <- p1 / p2 / p3 / p4 / p5 / p6 / p7 + plot_layout(heights = c(3,4,1,3,3,3,3))

    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/gene_analysis/top_amp_oncogenes.png")
    ggsave(filename, p, w = 13, h = 8)

    # case examples of top 10 amplified oncogenes
    toPlot <- oncogene_df[oncogene_df$oncogene %in% onco_order[1:10],]
    toPlot$oncogene <- factor(toPlot$oncogene, levels = onco_order[1:10])
    toPlot$Classification <- factor(toPlot$Classification, levels = names(amplicon_class_pal))

    p <- ggplot(toPlot, aes(x = Classification, y = log2(Feature_median_copy_number), fill = Classification)) +
        geom_violin() + geom_boxplot(width = 0.2) +
        facet_wrap(~oncogene, ncol=5, scales = "free_y") +
        scale_fill_manual(values = amplicon_class_pal) +
        scale_x_discrete(labels = c("BFB", "Complex\nnon-cyclic", "Linear", "ecDNA")) +
        theme_bw() +
        theme(legend.position = "none", axis.title.x = element_blank())
    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/gene_analysis/top_amp_oncogene_cases.png")
    ggsave(filename, p, width = 11, height = 5)
}

#' Plot top oncogenes by amplicon class
#' 
plot_top_amp_oncogenes_class <- function(oncogene_df, class) {

    toPlot <- oncogene_df[oncogene_df$Classification == class,]
    count_gene <- data.frame(table(toPlot$oncogene))
    count_gene <- count_gene[order(count_gene$Freq, decreasing = TRUE),]
    top_genes <- as.character(count_gene$Var1[1:30])

    toPlot <- toPlot[toPlot$oncogene %in% top_genes,]
    toPlot$oncogene <- factor(toPlot$oncogene, levels = top_genes)
    count_gene$Var1 <- factor(count_gene$Var1, levels = top_genes)


    # plot count of oncogenes
    p1 <- ggplot(count_gene[1:30,], aes(x = Var1, y = Freq)) +
        geom_bar(stat = "identity") +
        geom_text(
                aes(label = Freq),
                vjust = -0.3,
                color = "black",
                size = 3
            ) +
        scale_y_continuous(limits = c(0, max(count_gene$Freq)+4)) +
        theme_bw() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank()
        ) +
        labs(y = "No. Amplicons\nw Oncogene")

    p2 <- ggplot(toPlot, aes(x = oncogene, y = log2(Feature_median_copy_number))) +
        geom_violin(fill = amplicon_class_pal[class]) +
        geom_boxplot(width = 0.2, fill = amplicon_class_pal[class]) +
        theme_bw() + 
        theme(
            axis.text.x = element_text(angle = 25, hjust = 1),
            axis.title.x = element_blank()
        ) +
        labs(y = "log2(Median CN)")

    p <- p1 / p2 + plot_layout(height = c(1, 2))
    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/gene_analysis/amp_oncogenes_", class, ".png")
    ggsave(filename, p, width = 11, height = 4)
}

#' Plot top genes across amplicon classes and CN
#' displays if gene is oncogene
#' 
plot_top_amp_genes <- function(top_gene, all_genes_df) {

    gene_order <- top_gene$Var1
    top_gene$Var1 <- factor(top_gene$Var1, levels = gene_order)
    top_gene$oncogene <- factor(top_gene$oncogene, levels = c("Oncogene", "Not"))

    # plot count of oncogenes
    p1 <- ggplot(top_gene, aes(x = Var1, y = Freq, fill = oncogene)) +
        geom_bar(stat = "identity", color = "black") +
        geom_text(
                aes(label = Freq),
                vjust = -0.3,
                color = "black",
                size = 3
            ) +
        scale_y_continuous(limits = c(0, 1000)) +
        scale_fill_manual(values = binary_pal) +
        theme_bw() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank()
        ) +
        labs(y = "No. Amplicons\nw Gene")

    # plot proportion of amplicon class with oncogene
    toPlot <- all_genes_df[all_genes_df$gene %in% gene_order,]
    toPlot <- toPlot %>%
        group_by(gene, Classification) %>%
        summarise(total = n(), .groups = "drop") %>%
        group_by(gene) %>%
        mutate(prop = total / sum(total) * 100)
    toPlot$gene <- factor(toPlot$gene, levels = gene_order)
    toPlot$Classification <- factor(toPlot$Classification, levels = names(amplicon_class_pal))
    p2 <- ggplot(toPlot, aes(x = gene, y = prop, fill = Classification)) +
        geom_col(color = "black") +
        geom_text(
            aes(label = total),
            position = position_stack(vjust = 0.5),
            color = "black",
            size = 3
        ) +
        scale_fill_manual(values = amplicon_class_pal) + 
        theme_bw() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank()
        ) +
        labs(x = "Gene", y = "Proportion of\namplicons (%)", fill = "Amplicon\nClass")

    # helper function to plot distribution of feature median CN
    helper_plot_fmcn_class <- function(df, class) {

        toPlot <- df[df$gene %in% gene_order,]
        toPlot <- toPlot[toPlot$Classification == class,]
        toPlot$gene <- factor(toPlot$gene, levels = gene_order)

        p <- ggplot(toPlot, aes(x = gene, y = log2(Feature_median_copy_number))) +
            geom_boxplot(fill = amplicon_class_pal[class]) +
            scale_x_discrete(drop = FALSE) +
            #scale_y_continuous(limits = c(2, 7.5)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
            labs(y = "log2(Median\nCopy Number)", x = "Gene")
        if (class != "ecDNA") {
            p <- p + theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank()
            )
        }
        return(p)
    }

    p3 <- helper_plot_fmcn_class(all_genes_df, "BFB")
    p4 <- helper_plot_fmcn_class(all_genes_df, "Complex-non-cyclic")
    p5 <- helper_plot_fmcn_class(all_genes_df, "Linear")
    p6 <- helper_plot_fmcn_class(all_genes_df, "ecDNA")

    p <- p1 / p2 / p3 / p4 / p5 / p6 + plot_layout(heights = c(3,4,3,3,3,3))

    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/gene_analysis/top_amp_all_genes.png")
    ggsave(filename, p, w = 15, h = 8)
}

#' Plot top oncogenes by amplicon class
#' 
plot_top_amp_all_genes_class <- function(all_genes_df, class) {

    toPlot <- all_genes_df[all_genes_df$Classification == class,]
    count_gene <- data.frame(table(toPlot$gene))
    count_gene <- count_gene[order(count_gene$Freq, decreasing = TRUE),]
    top_genes <- as.character(count_gene$Var1[1:30])

    toPlot <- toPlot[toPlot$gene %in% top_genes,]
    toPlot$gene <- factor(toPlot$gene, levels = top_genes)
    count_gene$Var1 <- factor(count_gene$Var1, levels = top_genes)
    count_gene$oncogene <- toPlot$oncogene[match(count_gene$Var1, toPlot$gene)]
    count_gene$oncogene <- factor(count_gene$oncogene, levels = c("Oncogene", "Not"))

    # plot count of oncogenes
    p1 <- ggplot(count_gene[1:30,], aes(x = Var1, y = Freq, fill = oncogene)) +
        geom_bar(stat = "identity", color = "black") +
        geom_text(
                aes(label = Freq),
                vjust = -0.3,
                color = "black",
                size = 3
            ) +
        scale_fill_manual(values = binary_pal) +
        scale_y_continuous(limits = c(0, max(count_gene$Freq)+25)) +
        theme_bw() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank()
        ) +
        labs(y = "No. Amplicons\nw Gene")

    p2 <- ggplot(toPlot, aes(x = gene, y = log2(Feature_median_copy_number))) +
        geom_violin(fill = amplicon_class_pal[class]) +
        geom_boxplot(width = 0.2, fill = amplicon_class_pal[class]) +
        theme_bw() + 
        theme(
            axis.text.x = element_text(angle = 25, hjust = 1),
            axis.title.x = element_blank()
        ) +
        labs(y = "log2(Median CN)")

    p <- p1 / p2 + plot_layout(height = c(1, 2))
    filename <- paste0(RESULTS_DIR, "figures/AA_exploration/gene_analysis/amp_all_genes_", class, ".png")
    ggsave(filename, p, width = 11, height = 4)
}
