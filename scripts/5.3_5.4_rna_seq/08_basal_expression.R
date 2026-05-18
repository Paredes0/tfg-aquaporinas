#!/usr/bin/env Rscript
# =============================================================================
# 08_basal_expression.R — Basal Aquaporin Expression Analysis Across Tissues
#   Analyzes expression of ~145 aquaporins across 6 tissues:
#   Green Fruit, Red Fruit, Crown, Leaf (control), Roots (control), Aux Bud
#
# SCIENTIFIC NOTE — AuxBud_1 limitation:
#   AuxBud (Auxiliary Bud) has only N=1 replicate (SRR11806824).
#   With N=1, variance estimation is impossible, making statistical comparisons
#   unreliable. AuxBud is therefore:
#     - INCLUDED in TPM matrices and heatmaps as descriptive reference
#     - EXCLUDED from DESeq2 normalization (which requires ≥2 replicates for
#       reliable size factor estimation within a group)
#     - TPM values for AuxBud are calculated independently using raw counts
#     - Clearly flagged in all output tables and plots as "N=1, descriptive only"
#
# FIXES applied:
#   - Uses absolute paths derived from script location
#   - AuxBud handled as descriptive reference (N=1 limitation documented)
#   - Dynamic outlier exclusion
# =============================================================================

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(tidyr)
    library(dplyr)
})

# ---- Configuration (absolute paths) -----------------------------------------
PROJECT_DIR  <- getwd()
COUNTS_DIR   <- file.path(PROJECT_DIR, "counts")
DESIGN_DIR   <- file.path(PROJECT_DIR, "design")
RESULTS_DIR  <- file.path(PROJECT_DIR, "results", "basal_aquaporins")
AQP_TABLE    <- file.path(PROJECT_DIR, "tabla_Aquaporinas_traduccion.tabular")
GENE_LENGTHS <- file.path(COUNTS_DIR, "gene_lengths.tsv")

# Outlier handling
OUTLIER_HANDLING <- Sys.getenv("OUTLIER_HANDLING", unset = "include")
OUTLIER_SAMPLE   <- Sys.getenv("OUTLIER_SAMPLE",   unset = "RootsCtrl_2")

TPM_THRESHOLD <- 1  # TPM threshold for considering a gene "expressed"

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

message("# =============================================================")
message("# BASAL AQUAPORIN EXPRESSION ANALYSIS")
message("# =============================================================")

# ---- Load data ---------------------------------------------------------------
counts_all <- read.csv(file.path(COUNTS_DIR, "counts_basal.csv"),
                   row.names = 1, check.names = FALSE)
design_all <- read.csv(file.path(DESIGN_DIR, "design_basal.csv"))
gene_lengths <- read.delim(GENE_LENGTHS)
aqp_info <- read.delim(AQP_TABLE, sep = "\t", stringsAsFactors = FALSE)

# Restricción a las 121 funcionales (8 candidatas a reanotación excluidas;
# la reanotación no se aborda en el TFG). El filtro depende del esquema:
#   - si existe la columna explícita `needs_reannotation` → se usa directamente
#   - si no, se reconstruye desde `fuente_seq` (GFF3_FALLBACK / MAKER_GFF3 = TRUE)
# La cuantificación featureCounts y la normalización DESeq2 se mantienen sobre
# el genoma completo; sólo se restringe el subconjunto de acuaporinas.
if ("needs_reannotation" %in% colnames(aqp_info)) {
    aqp_info$needs_reannotation <- as.logical(aqp_info$needs_reannotation)
} else {
    aqp_info$needs_reannotation <- aqp_info$fuente_seq %in% c("GFF3_FALLBACK", "MAKER_GFF3")
}
n_pre_filter <- nrow(aqp_info)
aqp_info <- aqp_info[!aqp_info$needs_reannotation, ]
message(paste0("# Acuaporinas funcionales tras filtro: ", nrow(aqp_info),
               "/", n_pre_filter, " (excluidas ", n_pre_filter - nrow(aqp_info),
               " candidatas a reanotación)"))

# Dynamic outlier exclusion
if (OUTLIER_HANDLING == "exclude" && OUTLIER_SAMPLE %in% design_all$sample) {
    message(paste0("# EXCLUDING outlier: ", OUTLIER_SAMPLE))
    design_all <- design_all[design_all$sample != OUTLIER_SAMPLE, ]
    counts_all <- counts_all[, design_all$sample, drop = FALSE]
}

message(paste0("# Counts: ", nrow(counts_all), " genes x ", ncol(counts_all), " samples"))
message(paste0("# Design: ", nrow(design_all), " samples across ", length(unique(design_all$tissue)), " tissues"))
message(paste0("# Aquaporins: ", nrow(aqp_info), " genes"))

# ---- Separate AuxBud (N=1) from replicated tissues --------------------------
AUXBUD_SAMPLES <- design_all$sample[design_all$tissue == "aux_bud"]
has_auxbud <- length(AUXBUD_SAMPLES) > 0

if (has_auxbud) {
    message("# NOTE: AuxBud has N=1 replicate — treated as descriptive reference only")
    design_rep <- design_all[design_all$tissue != "aux_bud", ]
    counts_rep <- counts_all[, design_rep$sample, drop = FALSE]

    design_auxbud <- design_all[design_all$tissue == "aux_bud", ]
    counts_auxbud <- counts_all[, design_auxbud$sample, drop = FALSE]
} else {
    design_rep <- design_all
    counts_rep <- counts_all
}

# ---- DESeq2 normalization (replicated tissues only) --------------------------
design_rep$tissue <- factor(design_rep$tissue)
dds <- DESeqDataSetFromMatrix(
    countData = round(counts_rep),
    colData   = design_rep,
    design    = ~ tissue
)

# Filter lowly expressed genes, but ALWAYS KEEP all aquaporins
is_aqp <- rownames(dds) %in% aqp_info$gene_id
keep <- rowSums(counts(dds)) >= 10 | is_aqp
dds <- dds[keep, ]

# Estimate size factors for normalization
dds <- estimateSizeFactors(dds)
norm_counts_rep <- counts(dds, normalized = TRUE)

# VST for visualization (replicated tissues)
vsd <- vst(dds, blind = TRUE)

message(paste0("# Genes after filtering (replicated tissues): ", nrow(dds)))

# ---- Calculate TPM -----------------------------------------------------------
calculate_tpm <- function(counts_matrix, lengths_df) {
    # Match gene lengths to count matrix
    gene_len <- lengths_df$length[match(rownames(counts_matrix), lengths_df$gene_id)]
    # Replace NAs with median length
    gene_len[is.na(gene_len)] <- median(gene_len, na.rm = TRUE)

    # RPK = reads per kilobase
    rpk <- counts_matrix / (gene_len / 1000)
    # TPM = RPK / sum(RPK) * 1e6
    tpm <- sweep(rpk, 2, colSums(rpk) / 1e6, "/")
    return(tpm)
}

# TPM for replicated tissues
tpm_rep <- calculate_tpm(counts_rep[rownames(norm_counts_rep), ], gene_lengths)

# TPM for AuxBud (calculated independently — no DESeq2 normalization possible)
if (has_auxbud) {
    common_genes <- intersect(rownames(tpm_rep), rownames(counts_auxbud))
    tpm_auxbud <- calculate_tpm(counts_auxbud[common_genes, , drop = FALSE], gene_lengths)
    # Combine into full TPM matrix
    tpm_all <- cbind(tpm_rep[common_genes, ], tpm_auxbud[common_genes, , drop = FALSE])
    # Also build combined normalized counts (use raw for AuxBud since no size factor)
    norm_counts_all <- cbind(
        norm_counts_rep[common_genes, ],
        as.matrix(counts_auxbud[common_genes, , drop = FALSE])  # raw counts as proxy
    )
    design_combined <- rbind(design_rep, design_auxbud)
} else {
    tpm_all <- tpm_rep
    norm_counts_all <- norm_counts_rep
    design_combined <- design_rep
    common_genes <- rownames(norm_counts_rep)
}

# ---- Extract aquaporin data --------------------------------------------------
aqp_ids <- aqp_info$gene_id
aqp_in_counts <- aqp_ids[aqp_ids %in% common_genes]
aqp_missing <- aqp_ids[!aqp_ids %in% common_genes]

message(paste0("# Aquaporins in filtered counts: ", length(aqp_in_counts), "/", length(aqp_ids)))
if (length(aqp_missing) > 0) {
    message(paste0("# Missing/filtered aquaporins: ", paste(aqp_missing, collapse = ", ")))
}

# Extract aquaporin-specific matrices
aqp_norm <- norm_counts_all[aqp_in_counts, , drop = FALSE]
aqp_tpm  <- tpm_all[aqp_in_counts, , drop = FALSE]
aqp_vsd  <- assay(vsd)[intersect(aqp_in_counts, rownames(assay(vsd))), , drop = FALSE]

# ---- Merge with aquaporin metadata ------------------------------------------
aqp_meta <- aqp_info[aqp_info$gene_id %in% aqp_in_counts, ]
aqp_meta <- aqp_meta[match(aqp_in_counts, aqp_meta$gene_id), ]

# Create display names
aqp_meta$display_name <- paste0(aqp_meta$aqp_family_subfamily, " (", aqp_meta$gene_id, ")")
# `needs_reannotation` ya viene en aqp_info tras el filtro inicial; debería ser FALSE
# en todas las filas restantes. Se mantiene para compatibilidad con anotaciones aguas abajo.
if (!"needs_reannotation" %in% colnames(aqp_meta)) {
    aqp_meta$needs_reannotation <- FALSE
}

# ---- Save expression matrices ------------------------------------------------

# TPM matrix
aqp_tpm_df <- as.data.frame(aqp_tpm)
aqp_tpm_df$gene_id <- rownames(aqp_tpm_df)
aqp_tpm_df <- merge(aqp_meta[, c("gene_id", "aqp_family_subfamily", "subfamilia_phylo",
                                  "fuente_seq", "needs_reannotation")],
                     aqp_tpm_df, by = "gene_id")
write.csv(aqp_tpm_df, file.path(RESULTS_DIR, "basal_aquaporins_tpm.csv"), row.names = FALSE)

# Normalized counts
aqp_norm_df <- as.data.frame(aqp_norm)
aqp_norm_df$gene_id <- rownames(aqp_norm_df)
aqp_norm_df <- merge(aqp_meta[, c("gene_id", "aqp_family_subfamily", "subfamilia_phylo")],
                      aqp_norm_df, by = "gene_id")
write.csv(aqp_norm_df, file.path(RESULTS_DIR, "basal_aquaporins_normalized.csv"), row.names = FALSE)

# ---- Summary table: mean TPM per tissue per aquaporin ------------------------
tpm_long <- as.data.frame(aqp_tpm) %>%
    mutate(gene_id = rownames(aqp_tpm)) %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "tpm") %>%
    left_join(design_combined, by = "sample") %>%
    group_by(gene_id, tissue) %>%
    summarise(
        mean_tpm = mean(tpm),
        sd_tpm   = sd(tpm),
        n        = n(),
        .groups  = "drop"
    )

# Flag tissues with N=1 (no SD possible)
tpm_long <- tpm_long %>%
    mutate(
        sd_tpm = ifelse(n == 1, NA_real_, sd_tpm),
        note   = ifelse(n == 1, "N=1_descriptive_only", "")
    )

# Add aquaporin metadata
tpm_summary <- tpm_long %>%
    left_join(aqp_meta[, c("gene_id", "aqp_family_subfamily", "subfamilia_phylo",
                            "fuente_seq", "TMHs", "needs_reannotation")],
              by = "gene_id")

write.csv(tpm_summary, file.path(RESULTS_DIR, "basal_aquaporins_summary.csv"), row.names = FALSE)

# ---- Detection table: expressed (TPM > threshold) per tissue -----------------
detection <- tpm_long %>%
    mutate(expressed = mean_tpm > TPM_THRESHOLD) %>%
    select(gene_id, tissue, expressed) %>%
    pivot_wider(names_from = tissue, values_from = expressed)

detection <- merge(aqp_meta[, c("gene_id", "aqp_family_subfamily", "subfamilia_phylo")],
                   detection, by = "gene_id")
write.csv(detection, file.path(RESULTS_DIR, "basal_aquaporins_detection.csv"), row.names = FALSE)

message(paste0("# Expressed aquaporins (TPM > ", TPM_THRESHOLD, ") per tissue:"))
det_cols <- setdiff(colnames(detection), c("gene_id", "aqp_family_subfamily", "subfamilia_phylo"))
for (tis in det_cols) {
    n_expr <- sum(detection[[tis]], na.rm = TRUE)
    suffix <- ifelse(tis == "aux_bud", " (N=1, descriptive)", "")
    message(paste0("#   ", tis, ": ", n_expr, "/", nrow(detection), suffix))
}

# ==============================================================================
# VISUALIZATIONS
# ==============================================================================

# ---- Tissue color palette ----------------------------------------------------
tissue_colors <- c(
    "green_fruit" = "#228B22",
    "red_fruit"   = "#DC143C",
    "crown"       = "#DAA520",
    "leaf"        = "#32CD32",
    "roots"       = "#8B4513",
    "aux_bud"     = "#9370DB"
)

subfamily_colors <- c(
    "PIP" = "#E74C3C", "TIP" = "#3498DB", "NIP" = "#2ECC71",
    "SIP" = "#F39C12", "XIP" = "#9B59B6"
)

# ---- 1. PCA of aquaporin expression (replicated tissues only) ----------------
# Note: PCA uses VST data from replicated tissues. AuxBud excluded from PCA
# because N=1 cannot contribute meaningful variance information.
aqp_vsd_mat <- aqp_vsd
# Remove zero-variance genes before PCA
var_genes <- apply(aqp_vsd_mat, 1, var)
aqp_vsd_mat <- aqp_vsd_mat[var_genes > 1e-6, , drop = FALSE]

if (ncol(aqp_vsd_mat) >= 3 && nrow(aqp_vsd_mat) >= 2) {
    pca <- prcomp(t(aqp_vsd_mat), scale. = TRUE)
    pca_df <- data.frame(
        PC1 = pca$x[, 1], PC2 = pca$x[, 2],
        sample = rownames(pca$x)
    )
    pca_df <- merge(pca_df, design_rep, by = "sample")
    pct_var <- round(100 * summary(pca)$importance[2, 1:2])

    p_pca <- ggplot(pca_df, aes(PC1, PC2, color = tissue, label = sample)) +
        geom_point(size = 4) +
        geom_text(vjust = -0.8, size = 2.5) +
        xlab(paste0("PC1: ", pct_var[1], "% variance")) +
        ylab(paste0("PC2: ", pct_var[2], "% variance")) +
        ggtitle("PCA of Aquaporin Expression Across Tissues (N≥2 only)") +
        labs(caption = "Note: AuxBud (N=1) excluded from PCA") +
        theme_bw() +
        scale_color_manual(values = tissue_colors)
    ggsave(file.path(RESULTS_DIR, "pca_aquaporins.pdf"), p_pca, width = 10, height = 7)
}

# ---- 2. Sample-to-sample correlation (all samples, AuxBud flagged) ----------
cor_mat <- cor(aqp_norm, method = "pearson")
anno_col <- data.frame(tissue = design_combined$tissue, row.names = design_combined$sample)

pdf(file.path(RESULTS_DIR, "correlation_samples.pdf"), width = 10, height = 8)
pheatmap(
    cor_mat,
    annotation_col = anno_col,
    annotation_row = anno_col,
    annotation_colors = list(tissue = tissue_colors),
    color = colorRampPalette(brewer.pal(9, "Blues"))(100),
    main = "Sample-to-Sample Correlation (Aquaporins)\n[AuxBud_1: N=1, descriptive only]",
    display_numbers = TRUE,
    number_format = "%.2f",
    fontsize_number = 6
)
dev.off()

# ---- 3. Main heatmap: all expressed aquaporins across tissues ----------------
# Use mean TPM per tissue for a cleaner heatmap
tpm_wide <- tpm_long %>%
    select(gene_id, tissue, mean_tpm) %>%
    pivot_wider(names_from = tissue, values_from = mean_tpm) %>%
    as.data.frame()
rownames(tpm_wide) <- tpm_wide$gene_id
tpm_wide$gene_id <- NULL

# Do not filter out low expressed aquaporins, show all of them.
tpm_expressed <- tpm_wide

if (nrow(tpm_expressed) > 1) {
    # log2(TPM+1) transform
    mat_log <- log2(tpm_expressed + 1)

    # Row annotation
    meta_match <- aqp_meta[match(rownames(mat_log), aqp_meta$gene_id), ]
    rownames(mat_log) <- meta_match$display_name

    anno_row <- data.frame(
        subfamily = meta_match$subfamilia_phylo,
        reannotation = ifelse(meta_match$needs_reannotation, "needs_check", "OK"),
        row.names = meta_match$display_name
    )

    # Column annotation: flag AuxBud as N=1
    col_anno <- data.frame(
        replicate_status = ifelse(colnames(mat_log) == "aux_bud", "N=1", "N≥2"),
        row.names = colnames(mat_log)
    )

    anno_colors <- list(
        subfamily = subfamily_colors,
        reannotation = c("OK" = "grey90", "needs_check" = "orange"),
        replicate_status = c("N=1" = "red", "N≥2" = "grey90")
    )

    pdf(file.path(RESULTS_DIR, "heatmap_aquaporins_basal.pdf"),
        width = 12, height = max(14, nrow(mat_log) * 0.2))
    pheatmap(
        mat_log,
        annotation_row = anno_row,
        annotation_col = col_anno,
        annotation_colors = anno_colors,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = TRUE,
        fontsize_row = 5.5,
        main = "Aquaporin Basal Expression Across Tissues (log2(TPM+1))\n[AuxBud: N=1, descriptive reference]",
        color = colorRampPalette(c("white", "#FFFFCC", "#FEB24C", "#F03B20", "#BD0026"))(100),
        border_color = NA
    )
    dev.off()

    # Scaled version (z-score per gene)
    mat_scaled <- t(scale(t(as.matrix(mat_log))))
    mat_scaled[is.nan(mat_scaled)] <- 0

    pdf(file.path(RESULTS_DIR, "heatmap_aquaporins_basal_scaled.pdf"),
        width = 12, height = max(14, nrow(mat_scaled) * 0.2))
    pheatmap(
        mat_scaled,
        annotation_row = anno_row,
        annotation_col = col_anno,
        annotation_colors = anno_colors,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = TRUE,
        fontsize_row = 5.5,
        main = "Aquaporin Basal Expression Across Tissues (z-score)\n[AuxBud: N=1, descriptive reference]",
        color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
        border_color = NA
    )
    dev.off()
}

# ---- 4. Barplots per subfamily -----------------------------------------------
subfamily_summary <- tpm_summary %>%
    group_by(subfamilia_phylo, tissue) %>%
    summarise(
        total_tpm = sum(mean_tpm),
        mean_tpm  = mean(mean_tpm),
        n_genes   = n(),
        n_expressed = sum(mean_tpm > TPM_THRESHOLD),
        .groups = "drop"
    ) %>%
    mutate(tissue_label = ifelse(tissue == "aux_bud",
                                  paste0(tissue, "\n(N=1)"),
                                  tissue))

# Total expression per subfamily per tissue
p_bar_total <- ggplot(subfamily_summary,
                      aes(x = tissue_label, y = total_tpm, fill = tissue)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ subfamilia_phylo, scales = "free_y", ncol = 3) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = tissue_colors) +
    labs(title = "Total Aquaporin Expression per Subfamily per Tissue",
         subtitle = "AuxBud (N=1): shown as descriptive reference",
         x = "", y = "Sum TPM")
ggsave(file.path(RESULTS_DIR, "barplot_subfamily_total.pdf"),
       p_bar_total, width = 14, height = 10)

# Mean expression per subfamily per tissue
p_bar_mean <- ggplot(subfamily_summary,
                     aes(x = tissue_label, y = mean_tpm, fill = tissue)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ subfamilia_phylo, scales = "free_y", ncol = 3) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = tissue_colors) +
    labs(title = "Mean Aquaporin Expression per Subfamily per Tissue",
         subtitle = "AuxBud (N=1): shown as descriptive reference",
         x = "", y = "Mean TPM")
ggsave(file.path(RESULTS_DIR, "barplot_subfamily_mean.pdf"),
       p_bar_mean, width = 14, height = 10)

# ---- 5. Expression profile: dotplots per subfamily --------------------------
profile_data <- tpm_summary %>%
    filter(mean_tpm > 0.1) %>%
    mutate(tissue = factor(tissue, levels = c("green_fruit", "red_fruit", "crown",
                                               "leaf", "roots", "aux_bud")))

for (subfam in unique(profile_data$subfamilia_phylo)) {
    sub_data <- profile_data %>% filter(subfamilia_phylo == subfam)

    if (nrow(sub_data) == 0) next

    p_profile <- ggplot(sub_data,
                        aes(x = tissue, y = mean_tpm,
                            group = aqp_family_subfamily,
                            color = aqp_family_subfamily)) +
        geom_point(size = 3) +
        geom_line(linewidth = 0.8) +
        geom_errorbar(aes(ymin = pmax(mean_tpm - ifelse(is.na(sd_tpm), 0, sd_tpm), 0),
                          ymax = mean_tpm + ifelse(is.na(sd_tpm), 0, sd_tpm)),
                      width = 0.2, linewidth = 0.4) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "right") +
        scale_y_log10() +
        labs(title = paste0(subfam, " Aquaporin Expression Profiles"),
             subtitle = "Error bars: ±SD (absent for AuxBud, N=1)",
             x = "", y = "Mean TPM (log10 scale)",
             color = "Aquaporin")

    ggsave(file.path(RESULTS_DIR, paste0("profile_", subfam, ".pdf")),
           p_profile, width = 12, height = 7)
}

# ---- 6. Summary by number of expressed genes per tissue ----------------------
expressed_by_tissue <- tpm_summary %>%
    mutate(expressed = mean_tpm > TPM_THRESHOLD) %>%
    group_by(tissue, subfamilia_phylo) %>%
    summarise(
        n_expressed = sum(expressed),
        n_total = n(),
        .groups = "drop"
    ) %>%
    mutate(tissue_label = ifelse(tissue == "aux_bud",
                                  paste0(tissue, " (N=1)"),
                                  tissue))

p_expressed <- ggplot(expressed_by_tissue,
                      aes(x = tissue_label, y = n_expressed, fill = subfamilia_phylo)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = subfamily_colors) +
    labs(title = paste0("Number of Expressed Aquaporins per Tissue (TPM > ", TPM_THRESHOLD, ")"),
         subtitle = "AuxBud (N=1): descriptive reference, no variance estimate",
         x = "", y = "Number of aquaporins", fill = "Subfamily")
ggsave(file.path(RESULTS_DIR, "expressed_count_by_tissue.pdf"),
       p_expressed, width = 10, height = 6)

message("\n=== 08_basal_expression.R complete ===")
message(paste0("# Results saved to: ", RESULTS_DIR))
message("# IMPORTANT: AuxBud_1 (N=1) is included as descriptive reference only.")
message("# AuxBud_1 was excluded from DESeq2 normalization (requires N≥2).")
message("# AuxBud_1 TPM values were calculated independently from raw counts.")
