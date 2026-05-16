#!/usr/bin/env Rscript
# =============================================================================
# 14_homeolog_de_analysis.R
#   Performs DE analysis collapsing ONLY the homeolog aquaporin groups,
#   keeping the rest of the 100k+ genome features intact for proper dispersion
#   and size factor estimation in DESeq2.
#
#   Matches the exact output format and plot styles of 07_de_analysis.R
# =============================================================================

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(ggrepel)
    library(dplyr)
})

# ---- Configuration ----------------------------------------------------------
PROJECT_DIR <- getwd()
COUNTS_DIR  <- file.path(PROJECT_DIR, "counts")
DESIGN_DIR  <- file.path(PROJECT_DIR, "design")
RESULTS_DIR <- file.path(PROJECT_DIR, "results", "homeolog_de_analysis")
HG_FILE     <- file.path(PROJECT_DIR, "homeolog_groups.tsv")

OUTLIER_HANDLING <- Sys.getenv("OUTLIER_HANDLING", unset = "exclude")
OUTLIER_SAMPLE   <- Sys.getenv("OUTLIER_SAMPLE",   unset = "RootsCtrl_2")

FDR_CUTOFF <- 0.05
LFC_CUTOFF <- 1

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Color palettes exactly as in 07/08
sf_colors <- c(
    "PIP" = "#E74C3C", "TIP" = "#3498DB", "NIP" = "#2ECC71",
    "SIP" = "#F39C12", "XIP" = "#9B59B6", "Unknown" = "#999999"
)
cond_colors <- c("control" = "#4DAF4A", "stress" = "#E41A1C")

# Load data
counts_raw <- read.csv(file.path(COUNTS_DIR, "counts_all.csv"), row.names = 1, check.names = FALSE)
hg <- read.delim(HG_FILE, stringsAsFactors = FALSE)

# Build group annotation
hg_annot <- hg %>% group_by(homeolog_group) %>%
    summarise(family = first(family), sub_subfamily = first(sub_subfamily), .groups="drop")
rownames(hg_annot) <- hg_annot$homeolog_group

# Map gene -> ID (gene itself usually, or homeolog group if it's an aquaporin)
gene_to_grp <- setNames(rep(NA, nrow(counts_raw)), rownames(counts_raw))
names(gene_to_grp) <- rownames(counts_raw)
gene_to_grp[] <- rownames(counts_raw)

present_aqps <- hg$gene_id[hg$gene_id %in% rownames(counts_raw)]
gene_to_grp[present_aqps] <- hg$homeolog_group[match(present_aqps, hg$gene_id)]

message("# Collapsing matrix (summing homeologs)...")
collapsed_counts <- rowsum(as.matrix(counts_raw), group = gene_to_grp)
message(sprintf("# Collapsed matrix: %d rows (from %d original)", nrow(collapsed_counts), nrow(counts_raw)))

unique_hgs <- unique(hg$homeolog_group)

# ---- Helper function: run DE for a given comparison -------------------------
run_homeolog_de <- function(design_file, comparison_name) {
    message(paste0("\n# ========== ", comparison_name, " =========="))
    
    design <- read.csv(design_file, stringsAsFactors = FALSE)
    
    if (OUTLIER_HANDLING == "exclude" && OUTLIER_SAMPLE %in% design$sample) {
        message(paste0("# EXCLUDING outlier: ", OUTLIER_SAMPLE))
        design <- design[design$sample != OUTLIER_SAMPLE, ]
    }
    
    sub_counts <- collapsed_counts[, design$sample, drop=FALSE]
    design$condition <- factor(design$condition, levels = c("control", "stress"))
    
    # Run DESeq2
    dds <- DESeqDataSetFromMatrix(countData = round(sub_counts), colData = design, design = ~ condition)
    
    is_hg <- rownames(dds) %in% unique_hgs
    keep <- rowSums(counts(dds)) >= 10 | is_hg
    dds <- dds[keep, ]
    message(paste0("# Elements after filtering: ", nrow(dds)))
    
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", "stress", "control"))
    norm_counts <- counts(dds, normalized = TRUE)
    
    # Filter strictly to Homeolog Groups
    res_hg <- res[rownames(res) %in% unique_hgs, ]
    norm_hg <- norm_counts[rownames(norm_counts) %in% unique_hgs, , drop=FALSE]
    
    # Build Standardized Results DataFrame
    ctrl_samples   <- design$sample[design$condition == "control"]
    stress_samples <- design$sample[design$condition == "stress"]
    
    results_df <- data.frame(
        homeolog_group = rownames(res_hg),
        baseMean       = res_hg$baseMean,
        baseMeanA      = rowMeans(norm_hg[, ctrl_samples, drop = FALSE]),
        baseMeanB      = rowMeans(norm_hg[, stress_samples, drop = FALSE]),
        foldChange     = 2^res_hg$log2FoldChange,
        log2FoldChange = res_hg$log2FoldChange,
        PValue         = res_hg$pvalue,
        PAdj           = res_hg$padj,
        stringsAsFactors = FALSE
    )
    
    # Merge annot
    results_df <- merge(results_df, hg_annot, by="homeolog_group", all.x=TRUE)
    
    # Calculate FDR and sort
    results_df <- results_df[order(results_df$PAdj, na.last = TRUE), ]
    results_df$FDR <- results_df$PAdj
    results_df$falsePos <- seq_len(nrow(results_df)) * ifelse(is.na(results_df$FDR), NA, results_df$FDR)
    
    # Append normalized counts
    results_df <- cbind(results_df, as.data.frame(norm_hg[results_df$homeolog_group, , drop=FALSE]))
    
    out_csv <- file.path(RESULTS_DIR, paste0("de_homeologs_", comparison_name, ".csv"))
    write.csv(results_df, out_csv, row.names = FALSE)
    message(paste0("# Output CSV: ", out_csv))
    
    # Summary
    sig <- results_df[!is.na(results_df$FDR) & results_df$FDR < FDR_CUTOFF, ]
    sig_up   <- sig[sig$log2FoldChange > 0, ]
    sig_down <- sig[sig$log2FoldChange < 0, ]
    message(paste0("# Significant HG (FDR < ", FDR_CUTOFF, "): ", nrow(sig)))
    message(paste0("#   Up-regulated (stress > control):   ", nrow(sig_up)))
    message(paste0("#   Down-regulated (stress < control): ", nrow(sig_down)))

    # ---- 1. PCA plot (using ALL background genes for accurate scaling) --------
    vsd <- vst(dds, blind = TRUE)
    pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    pct_var <- round(100 * attr(pca_data, "percentVar"))

    p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
        geom_point(size = 4) +
        geom_text(vjust = -0.8, size = 3) +
        xlab(paste0("PC1: ", pct_var[1], "% variance")) +
        ylab(paste0("PC2: ", pct_var[2], "% variance")) +
        ggtitle(paste0("PCA — HG DE ", comparison_name)) +
        theme_bw() +
        scale_color_manual(values = cond_colors)
    ggsave(file.path(RESULTS_DIR, paste0("pca_", comparison_name, ".pdf")), p_pca, width = 8, height = 6)

    # ---- 2. MA plot ----------------------------------------------------------
    pdf(file.path(RESULTS_DIR, paste0("ma_", comparison_name, ".pdf")), width = 8, height = 6)
    plotMA(res, main = paste0("MA plot — ", comparison_name), ylim = c(-5, 5))
    dev.off()

    # ---- 3. Custom Volcano Plot (matching 07_de_analysis.R exactly with bg cloud) ---
    volcano_df <- as.data.frame(res)
    volcano_df$feature_id <- rownames(volcano_df)
    volcano_df$neg_log10_padj <- -log10(volcano_df$padj)
    volcano_df$neg_log10_padj[is.infinite(volcano_df$neg_log10_padj)] <- 
        max(volcano_df$neg_log10_padj[is.finite(volcano_df$neg_log10_padj)], na.rm = TRUE) + 5
    
    # Identify which are our Homeolog Groups
    volcano_df$is_hg <- volcano_df$feature_id %in% unique_hgs
    
    # Significance categories for non-hg background coloring
    volcano_df$signif_cat <- "NS"
    volcano_df$signif_cat[!is.na(volcano_df$padj) & volcano_df$padj < FDR_CUTOFF &
                           abs(volcano_df$log2FoldChange) >= LFC_CUTOFF] <- "Sig"
                           
    # Merge annot for foreground
    volcano_df <- merge(volcano_df, hg_annot, by.x="feature_id", by.y="homeolog_group", all.x=TRUE)
    
    # For labeling: only significant and only Homeolog Groups
    volcano_df$label <- NA_character_
    is_sig_hg <- volcano_df$is_hg & !is.na(volcano_df$padj) & volcano_df$padj < FDR_CUTOFF
    volcano_df$label[is_sig_hg] <- paste0(volcano_df$feature_id[is_sig_hg], " (", volcano_df$sub_subfamily[is_sig_hg], ")")
    
    bg_df <- volcano_df[!volcano_df$is_hg, ]
    fg_df <- volcano_df[volcano_df$is_hg, ]
    
    p_volc <- ggplot() +
        # Background: all normal genes (non-aquaporin)
        geom_point(data = bg_df,
            aes(x = log2FoldChange, y = neg_log10_padj),
            color = ifelse(bg_df$signif_cat == "Sig", "#AAAAAA", "#DDDDDD"),
            size = 0.4, alpha = 0.5) +
        # Foreground: Homeolog Groups marked by Family
        geom_point(data = fg_df,
            aes(x = log2FoldChange, y = neg_log10_padj, color = family),
            size = 3, alpha = 0.9) +
        # Labels for significant HG
        geom_text_repel(data = fg_df[!is.na(fg_df$label), ],
            aes(x = log2FoldChange, y = neg_log10_padj, label=label, color=family),
            size = 3.2, fontface = "bold", box.padding = 0.5, point.padding = 0.3, show.legend=FALSE) +
        # Threshold lines
        geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", color = "grey40", linewidth = 0.4) +
        geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", color = "grey40", linewidth = 0.4) +
        scale_color_manual(values = sf_colors, name = "AQ Family") +
        labs(
            title = paste0("Volcano — HG DE ", comparison_name),
            subtitle = "Stress vs Control · Homeolog Groups highlighted over full transcriptome",
            x = expression(log[2]~Fold~Change),
            y = expression(-log[10]~adjusted~p~value),
            caption = paste0("Labels: Significant Homeolog Groups | ",
                             nrow(fg_df), " HGs total | ",
                             "Dashed lines: FDR=", FDR_CUTOFF, ", |log2FC|=", LFC_CUTOFF)
        ) + theme_bw(base_size=12) +
        theme(
            plot.title = element_text(face="bold", size=14),
            plot.subtitle = element_text(color="grey40"),
            legend.position = "right"
        )
    
    ggsave(file.path(RESULTS_DIR, paste0("volcano_homeologs_", comparison_name, ".pdf")), p_volc, width = 12, height = 9)

    # ---- 4. Heatmap of ALL Homeolog Groups ----------------------------------
    mat <- norm_hg
    # Exclude rows with 0 variance
    var_genes <- apply(mat, 1, var)
    mat <- mat[var_genes > 1e-6, , drop=FALSE]
    
    if (nrow(mat) >= 2) {
        mat_scaled <- t(scale(t(log2(mat + 1))))
        mat_scaled[is.nan(mat_scaled)] <- 0
        
        # Add family label to rownames
        rn_mapped <- merge(data.frame(homeolog_group=rownames(mat_scaled)), hg_annot, by="homeolog_group", sort=FALSE)
        rownames(mat_scaled) <- paste0(rn_mapped$homeolog_group, " (", rn_mapped$sub_subfamily, ")")
        
        anno_row <- data.frame(family = rn_mapped$family, row.names=rownames(mat_scaled))
        anno_col <- data.frame(condition = design$condition, row.names = design$sample)

        pdf(file.path(RESULTS_DIR, paste0("heatmap_homeologs_", comparison_name, ".pdf")),
            width = 10, height = max(6, nrow(mat_scaled)*0.2))
        pheatmap(
            mat_scaled,
            annotation_col = anno_col,
            annotation_row = anno_row,
            annotation_colors = list(condition = cond_colors, family = sf_colors),
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = TRUE,
            fontsize_row = 6,
            main = paste0("Homeolog Group Expression — ", comparison_name, "\nRow-scaled log2(norm_counts + 1)"),
            color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
        )
        dev.off()
    }
}

# Run for both
run_homeolog_de(file.path(DESIGN_DIR, "design_leaf_de.csv"), "leaf")
run_homeolog_de(file.path(DESIGN_DIR, "design_roots_de.csv"), "roots")

message("\n# Done! Standardized DE results in ", RESULTS_DIR)
