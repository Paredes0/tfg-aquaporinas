#!/usr/bin/env Rscript
# =============================================================================
# 07_de_analysis.R — Differential Expression with DESeq2
#   Performs DE analysis for Leaf (control vs stress) and Roots (control vs stress)
#   Generates standardized output following Biostar Handbook conventions
#
# FIXES applied:
#   - Uses absolute paths derived from script location (not relative)
#   - Dynamic outlier exclusion based on environment variable
#   - Robust error handling
# =============================================================================

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
})

# ---- Configuration (absolute paths from script location) --------------------
# Determine project root from script location
script_dir <- tryCatch(
    dirname(normalizePath(sys.frame(1)$ofile)),
    error = function(e) getwd()
)
# If running via Rscript from PROJECT_DIR, getwd() is already correct
PROJECT_DIR <- getwd()

COUNTS_DIR   <- file.path(PROJECT_DIR, "counts")
DESIGN_DIR   <- file.path(PROJECT_DIR, "design")
RESULTS_DIR  <- file.path(PROJECT_DIR, "results")
AQP_TABLE    <- file.path(PROJECT_DIR, "tabla_Aquaporinas_traduccion.tabular")

# Read outlier handling from environment (set by run_pipeline.sh)
OUTLIER_HANDLING <- Sys.getenv("OUTLIER_HANDLING", unset = "include")
OUTLIER_SAMPLE   <- Sys.getenv("OUTLIER_SAMPLE",   unset = "RootsCtrl_2")

# FDR threshold for significance
FDR_CUTOFF <- 0.05
LFC_CUTOFF <- 1  # log2FC threshold for volcano plot coloring

# ---- Helper function: run DE for a given comparison -------------------------
run_de_analysis <- function(counts_file, design_file, output_dir, comparison_name, aqp_info) {

    message(paste0("# ========== ", comparison_name, " =========="))

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # Load counts
    counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
    message(paste0("# Counts: ", nrow(counts), " genes x ", ncol(counts), " samples"))

    # Load design
    design <- read.csv(design_file)
    message(paste0("# Design: ", nrow(design), " samples"))

    # Dynamic outlier exclusion
    if (OUTLIER_HANDLING == "exclude" && OUTLIER_SAMPLE %in% design$sample) {
        message(paste0("# EXCLUDING outlier: ", OUTLIER_SAMPLE))
        design <- design[design$sample != OUTLIER_SAMPLE, ]
        counts <- counts[, design$sample, drop = FALSE]
    }

    # Ensure sample order matches
    counts <- counts[, design$sample, drop = FALSE]

    # Set condition as factor with control as reference
    design$condition <- factor(design$condition, levels = c("control", "stress"))

    # Create DESeq dataset
    dds <- DESeqDataSetFromMatrix(
        countData = round(counts),  # ensure integers
        colData   = design,
        design    = ~ condition
    )

    # Pre-filter: remove genes with very low counts, EXCEPT aquaporins
    is_aqp <- rownames(dds) %in% aqp_info$gene_id
    keep <- rowSums(counts(dds)) >= 10 | is_aqp
    dds <- dds[keep, ]
    message(paste0("# Genes after filtering: ", nrow(dds)))

    # Run DESeq2
    dds <- DESeq(dds)

    # Extract results (stress vs control)
    res <- results(dds, contrast = c("condition", "stress", "control"))

    # Get normalized counts
    norm_counts <- counts(dds, normalized = TRUE)

    # ---- Build standardized results table (Biostar Handbook format) ----------
    ctrl_samples  <- design$sample[design$condition == "control"]
    stress_samples <- design$sample[design$condition == "stress"]

    results_df <- data.frame(
        name           = rownames(res),
        baseMean       = res$baseMean,
        baseMeanA      = rowMeans(norm_counts[, ctrl_samples, drop = FALSE]),
        baseMeanB      = rowMeans(norm_counts[, stress_samples, drop = FALSE]),
        foldChange     = 2^res$log2FoldChange,
        log2FoldChange = res$log2FoldChange,
        PValue         = res$pvalue,
        PAdj           = res$padj,
        stringsAsFactors = FALSE
    )

    # Add FDR and falsePos columns (ordered by PAdj)
    results_df <- results_df[order(results_df$PAdj, na.last = TRUE), ]
    results_df$FDR <- results_df$PAdj  # In DESeq2, padj is BH-adjusted = FDR
    results_df$falsePos <- seq_len(nrow(results_df)) * ifelse(
        is.na(results_df$FDR), NA, results_df$FDR
    )

    # Append normalized counts per sample
    norm_df <- as.data.frame(norm_counts[results_df$name, , drop = FALSE])
    results_df <- cbind(results_df, norm_df)

    # Save results
    results_file <- file.path(output_dir, paste0("results_", comparison_name, ".csv"))
    write.csv(results_df, results_file, row.names = FALSE)
    message(paste0("# Output: ", results_file))

    # ---- Summary statistics --------------------------------------------------
    sig <- results_df[!is.na(results_df$FDR) & results_df$FDR < FDR_CUTOFF, ]
    sig_up   <- sig[sig$log2FoldChange > 0, ]
    sig_down <- sig[sig$log2FoldChange < 0, ]
    message(paste0("# Significant genes (FDR < ", FDR_CUTOFF, "): ", nrow(sig)))
    message(paste0("#   Up-regulated (stress > control):   ", nrow(sig_up)))
    message(paste0("#   Down-regulated (stress < control): ", nrow(sig_down)))

    # ---- Plots ---------------------------------------------------------------

    # 1. PCA plot
    vsd <- vst(dds, blind = TRUE)
    pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    pct_var <- round(100 * attr(pca_data, "percentVar"))

    p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
        geom_point(size = 4) +
        geom_text(vjust = -0.8, size = 3) +
        xlab(paste0("PC1: ", pct_var[1], "% variance")) +
        ylab(paste0("PC2: ", pct_var[2], "% variance")) +
        ggtitle(paste0("PCA — ", comparison_name)) +
        theme_bw() +
        scale_color_manual(values = c("control" = "#4DAF4A", "stress" = "#E41A1C"))
    ggsave(file.path(output_dir, paste0("pca_", comparison_name, ".pdf")),
           p_pca, width = 8, height = 6)

    # 2. Volcano plot — only label aquaporin genes (custom ggplot2)
    suppressPackageStartupMessages(library(ggrepel))

    volcano_df <- as.data.frame(res)
    volcano_df$gene_id <- rownames(volcano_df)
    volcano_df$neg_log10_padj <- -log10(volcano_df$padj)
    # Cap infinite values
    volcano_df$neg_log10_padj[is.infinite(volcano_df$neg_log10_padj)] <- 
        max(volcano_df$neg_log10_padj[is.finite(volcano_df$neg_log10_padj)], na.rm = TRUE) + 5

    # Merge with aquaporin info
    volcano_df <- merge(volcano_df,
        aqp_info[, c("gene_id", "aqp_family_subfamily", "subfamilia_phylo")],
        by = "gene_id", all.x = TRUE
    )
    volcano_df$is_aqp <- !is.na(volcano_df$subfamilia_phylo)

    # Significance categories for non-aqp background coloring
    volcano_df$signif_cat <- "NS"
    volcano_df$signif_cat[!is.na(volcano_df$padj) & volcano_df$padj < FDR_CUTOFF &
                           abs(volcano_df$log2FoldChange) >= LFC_CUTOFF] <- "Sig"

    # Get top 20 most significant aquaporins for labeling
    aqp_df <- volcano_df[volcano_df$is_aqp, ]
    aqp_df <- aqp_df[order(aqp_df$padj, na.last = TRUE), ]
    top20_ids <- head(aqp_df$gene_id, 20)
    volcano_df$label <- ifelse(
        volcano_df$gene_id %in% top20_ids,
        volcano_df$aqp_family_subfamily,
        NA_character_
    )

    # Subfamily color palette
    sf_colors <- c(
        "PIP" = "#E74C3C", "TIP" = "#3498DB", "NIP" = "#2ECC71",
        "SIP" = "#F39C12", "XIP" = "#9B59B6"
    )

    # Split data: background (non-aqp) and foreground (aqp)
    bg_df <- volcano_df[!volcano_df$is_aqp, ]
    fg_df <- volcano_df[volcano_df$is_aqp, ]

    p_volcano <- ggplot() +
        # Background: all non-aquaporin genes as small grey dots
        geom_point(data = bg_df,
            aes(x = log2FoldChange, y = neg_log10_padj),
            color = ifelse(bg_df$signif_cat == "Sig", "#AAAAAA", "#DDDDDD"),
            size = 0.4, alpha = 0.5) +
        # Foreground: aquaporin genes colored by subfamily
        geom_point(data = fg_df,
            aes(x = log2FoldChange, y = neg_log10_padj, color = subfamilia_phylo),
            size = 3, alpha = 0.9) +
        # Labels for top 20 aquaporins
        geom_text_repel(data = fg_df[!is.na(fg_df$label), ],
            aes(x = log2FoldChange, y = neg_log10_padj, label = label, color = subfamilia_phylo),
            size = 3.2, fontface = "bold",
            box.padding = 0.5, point.padding = 0.3,
            segment.color = "grey50", segment.size = 0.3,
            max.overlaps = 30, show.legend = FALSE) +
        # Threshold lines
        geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", color = "grey40", linewidth = 0.4) +
        geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", color = "grey40", linewidth = 0.4) +
        # Subfamily colors
        scale_color_manual(values = sf_colors, name = "Subfamily") +
        # Labels
        labs(
            title = paste0("Volcano — ", comparison_name),
            subtitle = "Stress vs Control · Aquaporins highlighted",
            x = expression(log[2]~Fold~Change),
            y = expression(-log[10]~adjusted~p~value),
            caption = paste0("Top 20 aquaporins labeled | ",
                             nrow(fg_df), " aquaporins total | ",
                             "Dashed lines: FDR=", FDR_CUTOFF, ", |log2FC|=", LFC_CUTOFF)
        ) +
        theme_bw(base_size = 12) +
        theme(
            plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(color = "grey40"),
            legend.position = "right"
        )

    ggsave(file.path(output_dir, paste0("volcano_", comparison_name, ".pdf")),
           p_volcano, width = 12, height = 9)

    # 3. MA plot
    pdf(file.path(output_dir, paste0("ma_", comparison_name, ".pdf")), width = 8, height = 6)
    plotMA(res, main = paste0("MA plot — ", comparison_name), ylim = c(-5, 5))
    dev.off()

    # 4. Heatmap of top DE genes
    if (nrow(sig) > 1) {
        top_n <- min(50, nrow(sig))
        top_genes <- head(sig$name, top_n)
        mat <- norm_counts[top_genes, , drop = FALSE]
        mat_scaled <- t(scale(t(log2(mat + 1))))

        anno_col <- data.frame(
            condition = design$condition,
            row.names = design$sample
        )
        anno_colors <- list(condition = c("control" = "#4DAF4A", "stress" = "#E41A1C"))

        pdf(file.path(output_dir, paste0("heatmap_top_de_", comparison_name, ".pdf")),
            width = 8, height = max(6, top_n * 0.2))
        pheatmap(
            mat_scaled,
            annotation_col = anno_col,
            annotation_colors = anno_colors,
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = TRUE,
            fontsize_row = 6,
            main = paste0("Top ", top_n, " DE genes — ", comparison_name),
            color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
        )
        dev.off()
    }

    # ---- Return results for aquaporin filtering ------------------------------
    return(list(results = results_df, norm_counts = norm_counts, dds = dds))
}

# ---- Helper: annotate aquaporins in DE results -------------------------------
annotate_aquaporins <- function(de_results, aqp_info, output_dir, comparison_name, norm_counts, design_file) {

    design <- read.csv(design_file)

    # Dynamic outlier exclusion
    if (OUTLIER_HANDLING == "exclude" && OUTLIER_SAMPLE %in% design$sample) {
        design <- design[design$sample != OUTLIER_SAMPLE, ]
    }

    # Filter DE results for aquaporin genes
    aqp_de <- de_results[de_results$name %in% aqp_info$gene_id, ]
    aqp_de <- merge(aqp_de, aqp_info, by.x = "name", by.y = "gene_id", all.x = TRUE)

    # Sort by FDR
    aqp_de <- aqp_de[order(aqp_de$FDR, na.last = TRUE), ]

    # Flag partial/mis-annotated
    aqp_de$needs_reannotation <- aqp_de$fuente_seq %in% c("GFF3_FALLBACK", "MAKER_GFF3")

    # Save
    out_file <- file.path(output_dir, paste0("de_aquaporins_", comparison_name, ".csv"))
    write.csv(aqp_de, out_file, row.names = FALSE)
    message(paste0("# Aquaporin DE results: ", out_file))

    # Summary
    aqp_sig <- aqp_de[!is.na(aqp_de$FDR) & aqp_de$FDR < FDR_CUTOFF, ]
    message(paste0("# Significant aquaporins (FDR < ", FDR_CUTOFF, "): ", nrow(aqp_sig)))

    if (nrow(aqp_sig) > 0) {
        message("# Significant aquaporins:")
        for (i in seq_len(nrow(aqp_sig))) {
            msg <- paste0("#   ", aqp_sig$aqp_family_subfamily[i],
                          " (", aqp_sig$name[i], ")",
                          " log2FC=", round(aqp_sig$log2FoldChange[i], 2),
                          " FDR=", formatC(aqp_sig$FDR[i], format = "e", digits = 2))
            if (aqp_sig$needs_reannotation[i]) msg <- paste0(msg, " [NEEDS REANNOTATION]")
            message(msg)
        }
    }

    # Heatmap of aquaporin expression: KEEP ALL aquaporins
    aqp_expressed <- aqp_de  # Remove filter so all 144 are printed
    if (nrow(aqp_expressed) > 1) {
        samples <- design$sample
        mat <- norm_counts[aqp_expressed$name, samples, drop = FALSE]
        mat_scaled <- t(scale(t(log2(mat + 1))))
        mat_scaled[is.nan(mat_scaled)] <- 0  # Fix NaN for 0-count genes

        # Row labels: use aquaporin names
        rn <- paste0(aqp_expressed$aqp_family_subfamily, " (", aqp_expressed$name, ")")
        rownames(mat_scaled) <- rn

        anno_col <- data.frame(
            condition = design$condition,
            row.names = design$sample
        )

        # Row annotation: subfamily
        anno_row <- data.frame(
            subfamily = aqp_expressed$subfamilia_phylo,
            row.names = rn
        )
        subfamily_colors <- c(
            "PIP" = "#E74C3C", "TIP" = "#3498DB", "NIP" = "#2ECC71",
            "SIP" = "#F39C12", "XIP" = "#9B59B6"
        )

        pdf(file.path(output_dir, paste0("heatmap_aquaporins_", comparison_name, ".pdf")),
            width = 10, height = max(14, nrow(aqp_expressed) * 0.2))
        pheatmap(
            mat_scaled,
            annotation_col = anno_col,
            annotation_row = anno_row,
            annotation_colors = list(
                condition = c("control" = "#4DAF4A", "stress" = "#E41A1C"),
                subfamily = subfamily_colors
            ),
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = TRUE,
            fontsize_row = 5.5,
            main = paste0("Aquaporin expression — ", comparison_name),
            color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
        )
        dev.off()
    }

    return(aqp_de)
}

# ==============================================================================
# MAIN
# ==============================================================================

# Load aquaporin metadata
aqp_info <- read.delim(AQP_TABLE, sep = "\t", stringsAsFactors = FALSE)
# Select relevant columns
aqp_cols <- c("gene_id", "aqp_family_subfamily", "subfamilia_phylo",
              "fuente_seq", "veredicto", "TMHs", "longitud_aa")
aqp_info <- aqp_info[, intersect(aqp_cols, colnames(aqp_info))]

# ---- Leaf DE analysis --------------------------------------------------------
message("\n#################################################################")
message("# LEAF DIFFERENTIAL EXPRESSION: Control vs Stress (6 days)")
message("#################################################################\n")

leaf_result <- run_de_analysis(
    counts_file     = file.path(COUNTS_DIR, "counts_leaf_de.csv"),
    design_file     = file.path(DESIGN_DIR, "design_leaf_de.csv"),
    output_dir      = file.path(RESULTS_DIR, "de_leaf"),
    comparison_name = "leaf",
    aqp_info        = aqp_info
)

leaf_aqp <- annotate_aquaporins(
    de_results      = leaf_result$results,
    aqp_info        = aqp_info,
    output_dir      = file.path(RESULTS_DIR, "de_leaf"),
    comparison_name = "leaf",
    norm_counts     = leaf_result$norm_counts,
    design_file     = file.path(DESIGN_DIR, "design_leaf_de.csv")
)

# ---- Roots DE analysis -------------------------------------------------------
message("\n#################################################################")
message("# ROOTS DIFFERENTIAL EXPRESSION: Control vs Stress (6 days)")
message("#################################################################\n")

roots_result <- run_de_analysis(
    counts_file     = file.path(COUNTS_DIR, "counts_roots_de.csv"),
    design_file     = file.path(DESIGN_DIR, "design_roots_de.csv"),
    output_dir      = file.path(RESULTS_DIR, "de_roots"),
    comparison_name = "roots",
    aqp_info        = aqp_info
)

roots_aqp <- annotate_aquaporins(
    de_results      = roots_result$results,
    aqp_info        = aqp_info,
    output_dir      = file.path(RESULTS_DIR, "de_roots"),
    comparison_name = "roots",
    norm_counts     = roots_result$norm_counts,
    design_file     = file.path(DESIGN_DIR, "design_roots_de.csv")
)

message("\n=== 07_de_analysis.R complete ===")
