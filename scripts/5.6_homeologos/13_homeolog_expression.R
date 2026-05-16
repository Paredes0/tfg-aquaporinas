#!/usr/bin/env Rscript
# =============================================================================
# 13_homeolog_expression.R — Three-Level Homeolog Expression Analysis
#
# Analyzes aquaporin expression at three resolution levels:
#   Level 1 (Individual):  All 144 genes separately (heatmap, PCA)
#   Level 2 (Collapsed):   Summed counts/TPMs per homeolog group (36 groups)
#   Level 3 (Dominance):   Subgenome dominance within each homeolog group
#
# SCIENTIFIC RATIONALE:
#   In polyploid organisms, homeologous genes (copies across subgenomes)
#   share high sequence similarity. Analyzing expression at individual level
#   captures sub-functionalization, while collapsed analysis reduces noise
#   from homeolog-to-homeolog variation. Dominance analysis identifies
#   whether one subgenome consistently contributes more expression,
#   revealing patterns of biased fractionation (common in Fragaria,
#   where the F. vesca-derived subgenome is often dominant).
#
# METHODOLOGY:
#   - Collapsed counts: sum(raw_counts) per homeolog group → DESeq2-compatible
#   - Collapsed TPMs: sum(TPMs) per homeolog group → visualization
#   - Dominance: TPM_subgenome / sum(TPM_group) per tissue
#   - Tandem duplicates within a subgenome are summed BEFORE dominance calc
#     to avoid inflating copy-number-driven "dominance"
#
# INPUT:
#   - homeolog_groups.tsv             (from 11_homeolog_grouping.py)
#   - counts/counts_basal.csv         (raw counts matrix, all genes)
#   - design/design_basal.csv         (experimental design)
#   - counts/gene_lengths.tsv         (gene lengths for TPM calculation)
#   - tabla_Aquaporinas_traduccion.tabular (aquaporin gene list)
#
# OUTPUT (in results/homeolog_analysis/):
#   - collapsed_counts.csv            (36 groups × N samples)
#   - collapsed_tpm.csv               (36 groups × N samples)
#   - dominance_by_tissue.csv         (subgenome proportions)
#   - dominance_heatmap.pdf           (visualization)
#   - collapsed_heatmap.pdf           (homeolog-level expression)
#   - subgenome_barplot.pdf           (stacked barplot of proportions)
#   - homeolog_pca.pdf                (PCA at collapsed level)
#   - level_comparison.pdf            (individual vs collapsed PCA)
#
# USAGE:
#   Rscript 13_homeolog_expression.R
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(tidyr)
    library(dplyr)
    library(reshape2)
})

# ---- Configuration ----------------------------------------------------------
PROJECT_DIR  <- getwd()
COUNTS_DIR   <- file.path(PROJECT_DIR, "counts")
DESIGN_DIR   <- file.path(PROJECT_DIR, "design")
RESULTS_DIR  <- file.path(PROJECT_DIR, "results", "homeolog_analysis")
HG_FILE      <- file.path(PROJECT_DIR, "homeolog_groups.tsv")
AQP_TABLE    <- file.path(PROJECT_DIR, "tabla_Aquaporinas_traduccion.tabular")
GENE_LENGTHS <- file.path(COUNTS_DIR, "gene_lengths.tsv")

OUTLIER_HANDLING <- Sys.getenv("OUTLIER_HANDLING", unset = "include")
OUTLIER_SAMPLE   <- Sys.getenv("OUTLIER_SAMPLE",   unset = "RootsCtrl_2")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Color palette (consistent with pipeline) --------------------------------
SUBGENOME_COLORS <- c(
    "A" = "#E41A1C",   # red
    "B" = "#377EB8",   # blue
    "C" = "#4DAF4A",   # green
    "D" = "#FF7F00"    # orange
)

FAMILY_COLORS <- c(
    "PIP" = "#1B9E77",
    "TIP" = "#D95F02",
    "NIP" = "#7570B3",
    "SIP" = "#E7298A",
    "XIP" = "#66A61E"
)

message("# =================================================================")
message("# HOMEOLOG EXPRESSION ANALYSIS — Three-Level Resolution")
message("# =================================================================")

# =============================================================================
# LOAD DATA
# =============================================================================

# Homeolog groups
hg <- read.delim(HG_FILE, stringsAsFactors = FALSE)
message(paste0("# Homeolog groups: ", length(unique(hg$homeolog_group)),
               " groups, ", nrow(hg), " genes"))

# Raw counts
counts_file <- file.path(COUNTS_DIR, "counts_basal.csv")
if (!file.exists(counts_file)) {
    message("# WARNING: counts_basal.csv not found. Trying counts_de.csv...")
    counts_file <- file.path(COUNTS_DIR, "counts_de.csv")
}
if (!file.exists(counts_file)) {
    stop("ERROR: No counts file found in ", COUNTS_DIR)
}
counts_raw <- read.csv(counts_file, row.names = 1, check.names = FALSE)

# Design
design_file <- file.path(DESIGN_DIR, "design_basal.csv")
if (!file.exists(design_file)) {
    design_file <- file.path(DESIGN_DIR, "design_de.csv")
}
design <- read.csv(design_file, stringsAsFactors = FALSE)

# Gene lengths
gene_lengths <- read.delim(GENE_LENGTHS)

# Outlier exclusion
if (OUTLIER_HANDLING == "exclude" && OUTLIER_SAMPLE %in% design$sample) {
    message(paste0("# EXCLUDING outlier: ", OUTLIER_SAMPLE))
    design <- design[design$sample != OUTLIER_SAMPLE, ]
    counts_raw <- counts_raw[, design$sample, drop = FALSE]
}

message(paste0("# Counts: ", nrow(counts_raw), " genes × ", ncol(counts_raw), " samples"))

# =============================================================================
# TPM CALCULATION FUNCTION
# =============================================================================

calculate_tpm <- function(counts_matrix, lengths_df) {
    gene_len <- lengths_df$length[match(rownames(counts_matrix), lengths_df$gene_id)]
    gene_len[is.na(gene_len)] <- median(gene_len, na.rm = TRUE)
    rpk <- counts_matrix / (gene_len / 1000)
    tpm <- sweep(rpk, 2, colSums(rpk) / 1e6, "/")
    return(tpm)
}

# Full TPM matrix (all genes)
tpm_all <- calculate_tpm(counts_raw, gene_lengths)

# =============================================================================
# LEVEL 1: Individual Gene Analysis (extract aquaporin subset)
# =============================================================================
message("\n# ── LEVEL 1: Individual Gene Analysis ──")

aqp_genes <- hg$gene_id
aqp_in_counts <- aqp_genes[aqp_genes %in% rownames(counts_raw)]
aqp_missing   <- aqp_genes[!aqp_genes %in% rownames(counts_raw)]

if (length(aqp_missing) > 0) {
    message(paste0("#   WARNING: ", length(aqp_missing),
                   " aquaporin genes not found in counts matrix"))
    message(paste0("#   Missing: ", paste(head(aqp_missing, 10), collapse = ", ")))
}
message(paste0("#   Using ", length(aqp_in_counts), " aquaporin genes"))

# Individual aquaporin counts & TPMs
counts_aqp <- counts_raw[aqp_in_counts, , drop = FALSE]
tpm_aqp    <- tpm_all[aqp_in_counts, , drop = FALSE]

# Annotate with group info
aqp_annotation <- hg[hg$gene_id %in% aqp_in_counts, ]
rownames(aqp_annotation) <- aqp_annotation$gene_id

# =============================================================================
# LEVEL 2: Collapsed Homeolog Group Analysis
# =============================================================================
message("\n# ── LEVEL 2: Collapsed Homeolog Group Analysis ──")

# Sum raw counts by homeolog group
collapse_by_group <- function(matrix_data, group_mapping) {
    # group_mapping: data.frame with gene_id and homeolog_group columns
    # Returns: matrix with homeolog_group as rownames
    groups <- unique(group_mapping$homeolog_group)
    result <- matrix(0, nrow = length(groups), ncol = ncol(matrix_data),
                     dimnames = list(groups, colnames(matrix_data)))

    for (g in groups) {
        member_genes <- group_mapping$gene_id[group_mapping$homeolog_group == g]
        present <- member_genes[member_genes %in% rownames(matrix_data)]
        if (length(present) == 1) {
            result[g, ] <- as.numeric(matrix_data[present, ])
        } else if (length(present) > 1) {
            result[g, ] <- colSums(matrix_data[present, , drop = FALSE])
        }
    }
    return(result)
}

# Collapse
collapsed_counts <- collapse_by_group(counts_aqp, hg[, c("gene_id", "homeolog_group")])
collapsed_tpm    <- collapse_by_group(tpm_aqp, hg[, c("gene_id", "homeolog_group")])

# Create annotation for collapsed groups
group_annotation <- hg %>%
    group_by(homeolog_group) %>%
    summarise(
        family = first(family),
        sub_subfamily = first(sub_subfamily),
        chr_number = first(chr_number),
        group_size = first(group_size),
        n_subgenomes = first(n_subgenomes),
        completeness = first(group_completeness),
        .groups = "drop"
    ) %>%
    as.data.frame()
rownames(group_annotation) <- group_annotation$homeolog_group

message(paste0("#   Collapsed ", length(aqp_in_counts), " genes → ",
               nrow(collapsed_counts), " homeolog groups"))
message(paste0("#   Complete quartets: ",
               sum(group_annotation$completeness == "complete")))

# Save collapsed matrices
write.csv(collapsed_counts, file.path(RESULTS_DIR, "collapsed_counts.csv"))
write.csv(collapsed_tpm,    file.path(RESULTS_DIR, "collapsed_tpm.csv"))
message("#   Saved collapsed_counts.csv and collapsed_tpm.csv")

# =============================================================================
# LEVEL 3: Subgenome Dominance Analysis
# =============================================================================
message("\n# ── LEVEL 3: Subgenome Dominance Analysis ──")

# For each homeolog group and each sample, calculate the proportion of
# expression contributed by each subgenome.
# IMPORTANT: Tandem duplicates within a subgenome are first summed together,
# then proportions are calculated at the subgenome level.

# Determine tissue for each sample
sample_tissue <- setNames(design$tissue, design$sample)

compute_dominance <- function(tpm_matrix, group_mapping) {
    # Returns: long-format data frame with columns:
    #   homeolog_group, sample, tissue, subgenome, tpm_subgenome, tpm_group, proportion

    results <- list()

    for (g in unique(group_mapping$homeolog_group)) {
        members <- group_mapping[group_mapping$homeolog_group == g, ]
        present <- members[members$gene_id %in% rownames(tpm_matrix), ]
        if (nrow(present) == 0) next

        for (samp in colnames(tpm_matrix)) {
            tissue <- sample_tissue[samp]

            # Sum TPM per subgenome (handles tandem duplicates)
            sg_tpm <- present %>%
                mutate(tpm = as.numeric(tpm_matrix[gene_id, samp])) %>%
                group_by(subgenome) %>%
                summarise(tpm_sg = sum(tpm, na.rm = TRUE), .groups = "drop")

            total_tpm <- sum(sg_tpm$tpm_sg)

            for (i in seq_len(nrow(sg_tpm))) {
                prop <- if (total_tpm > 0) sg_tpm$tpm_sg[i] / total_tpm else NA
                results[[length(results) + 1]] <- data.frame(
                    homeolog_group = g,
                    sample = samp,
                    tissue = tissue,
                    subgenome = sg_tpm$subgenome[i],
                    tpm_subgenome = sg_tpm$tpm_sg[i],
                    tpm_group = total_tpm,
                    proportion = prop,
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    do.call(rbind, results)
}

dominance_df <- compute_dominance(tpm_aqp, hg[, c("gene_id", "homeolog_group", "subgenome")])

# Add family/subfamily annotation
dominance_df <- dominance_df %>%
    left_join(group_annotation[, c("homeolog_group", "family", "sub_subfamily")],
              by = "homeolog_group")

# Average dominance by tissue (for visualization)
dominance_by_tissue <- dominance_df %>%
    filter(!is.na(proportion)) %>%
    group_by(homeolog_group, tissue, subgenome, family, sub_subfamily) %>%
    summarise(
        mean_proportion = mean(proportion),
        sd_proportion   = sd(proportion),
        mean_tpm        = mean(tpm_subgenome),
        .groups = "drop"
    )

# Overall dominance (across all tissues)
dominance_overall <- dominance_df %>%
    filter(!is.na(proportion)) %>%
    group_by(homeolog_group, subgenome, family, sub_subfamily) %>%
    summarise(
        mean_proportion = mean(proportion),
        sd_proportion   = sd(proportion),
        .groups = "drop"
    )

# Identify dominant subgenome per group
dominant_subgenome <- dominance_overall %>%
    group_by(homeolog_group) %>%
    slice_max(mean_proportion, n = 1) %>%
    ungroup() %>%
    select(homeolog_group, dominant_subgenome = subgenome,
           dominance_proportion = mean_proportion)

message(paste0("#   Dominance table: ", nrow(dominance_by_tissue), " rows"))
message(paste0("#   Dominant subgenomes: ",
               paste(table(dominant_subgenome$dominant_subgenome), collapse = ", ")))

# Save dominance tables
write.csv(dominance_by_tissue, file.path(RESULTS_DIR, "dominance_by_tissue.csv"),
          row.names = FALSE)
write.csv(dominance_overall, file.path(RESULTS_DIR, "dominance_overall.csv"),
          row.names = FALSE)
write.csv(dominant_subgenome, file.path(RESULTS_DIR, "dominant_subgenome.csv"),
          row.names = FALSE)
message("#   Saved dominance_by_tissue.csv, dominance_overall.csv, dominant_subgenome.csv")

# =============================================================================
# VISUALIZATIONS
# =============================================================================
message("\n# ── Generating Visualizations ──")

# ---- 1. Collapsed Heatmap (Level 2) -----------------------------------------
message("#   1/5 Collapsed heatmap...")

# Mean TPM per tissue for collapsed groups
tpm_mean_collapsed <- data.frame(row.names = rownames(collapsed_tpm))
for (tissue in unique(design$tissue)) {
    tissue_samples <- design$sample[design$tissue == tissue]
    tissue_samples <- tissue_samples[tissue_samples %in% colnames(collapsed_tpm)]
    if (length(tissue_samples) > 0) {
        tpm_mean_collapsed[[tissue]] <- rowMeans(
            collapsed_tpm[, tissue_samples, drop = FALSE])
    }
}

# Log2 transform for heatmap
heatmap_data <- log2(as.matrix(tpm_mean_collapsed) + 1)

# Row annotation
row_ann <- group_annotation[rownames(heatmap_data), c("family", "sub_subfamily"), drop = FALSE]

# Family color for annotation
ann_colors <- list(family = FAMILY_COLORS[unique(row_ann$family)])

tryCatch({
    pdf(file.path(RESULTS_DIR, "collapsed_heatmap.pdf"), width = 10, height = 14)
    pheatmap(heatmap_data,
             annotation_row = row_ann,
             annotation_colors = ann_colors,
             cluster_cols = FALSE,
             clustering_method = "ward.D2",
             color = colorRampPalette(c("white", "#FFF7EC", "#FEE8C8",
                                        "#FDD49E", "#FDBB84", "#FC8D59",
                                        "#EF6548", "#D7301F", "#990000"))(100),
             main = "Homeolog Group Expression (log2 TPM+1)\nCollapsed by summing homeologs",
             fontsize_row = 7,
             fontsize_col = 9,
             border_color = NA)
    dev.off()
    message("#     → collapsed_heatmap.pdf")
}, error = function(e) message(paste0("#     ERROR: ", e$message)))

# ---- 2. Subgenome Dominance Barplot (Level 3) --------------------------------
message("#   2/5 Subgenome dominance barplot...")

# Filter to complete/triplet groups with sufficient expression
expressed_groups <- dominance_overall %>%
    group_by(homeolog_group) %>%
    summarise(total_prop = sum(mean_proportion), .groups = "drop") %>%
    filter(total_prop > 0)

dom_plot_data <- dominance_overall %>%
    filter(homeolog_group %in% expressed_groups$homeolog_group) %>%
    mutate(homeolog_group = factor(homeolog_group,
        levels = group_annotation$homeolog_group[
            order(group_annotation$family, group_annotation$sub_subfamily)]))

tryCatch({
    p_dom <- ggplot(dom_plot_data,
                    aes(x = homeolog_group, y = mean_proportion,
                        fill = subgenome)) +
        geom_bar(stat = "identity", position = "stack", width = 0.8) +
        scale_fill_manual(values = SUBGENOME_COLORS, name = "Subgenome") +
        labs(title = "Subgenome Dominance per Homeolog Group",
             subtitle = "Proportion of total group expression contributed by each subgenome",
             x = "Homeolog Group",
             y = "Mean Proportion of Expression") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
              panel.grid.major.x = element_blank()) +
        geom_hline(yintercept = 0.25, linetype = "dashed", alpha = 0.3) +
        facet_grid(. ~ family, scales = "free_x", space = "free_x")

    ggsave(file.path(RESULTS_DIR, "subgenome_barplot.pdf"), p_dom,
           width = 16, height = 6)
    message("#     → subgenome_barplot.pdf")
}, error = function(e) message(paste0("#     ERROR: ", e$message)))

# ---- 3. Tissue-specific dominance heatmap ------------------------------------
message("#   3/5 Tissue-specific dominance heatmap...")

# Create a matrix: rows = homeolog_group:subgenome, cols = tissue
dom_matrix_data <- dominance_by_tissue %>%
    filter(homeolog_group %in% expressed_groups$homeolog_group) %>%
    mutate(row_id = paste0(homeolog_group, ":", subgenome)) %>%
    select(row_id, tissue, mean_proportion) %>%
    pivot_wider(names_from = tissue, values_from = mean_proportion) %>%
    as.data.frame()

rownames(dom_matrix_data) <- dom_matrix_data$row_id
dom_matrix_data$row_id <- NULL

# Row annotation for dominance heatmap
dom_row_ann <- data.frame(
    subgenome = sub(".*:", "", rownames(dom_matrix_data)),
    group = sub(":.*", "", rownames(dom_matrix_data)),
    stringsAsFactors = FALSE
)
dom_row_ann$family <- group_annotation[dom_row_ann$group, "family"]
rownames(dom_row_ann) <- rownames(dom_matrix_data)
dom_row_ann$group <- NULL

dom_ann_colors <- list(
    subgenome = SUBGENOME_COLORS,
    family = FAMILY_COLORS[unique(dom_row_ann$family)]
)

tryCatch({
    # Replace NA with 0 for clustering (NA = no expression → 0 dominance)
    dom_mat_clean <- as.matrix(dom_matrix_data)
    dom_mat_clean[is.na(dom_mat_clean)] <- 0
    dom_mat_clean[is.nan(dom_mat_clean)] <- 0
    dom_mat_clean[is.infinite(dom_mat_clean)] <- 0

    pdf(file.path(RESULTS_DIR, "dominance_heatmap.pdf"), width = 10, height = 20)
    pheatmap(dom_mat_clean,
             annotation_row = dom_row_ann,
             annotation_colors = dom_ann_colors,
             cluster_cols = FALSE,
             clustering_method = "ward.D2",
             color = colorRampPalette(c("#F7FBFF", "#DEEBF7", "#C6DBEF",
                                        "#9ECAE1", "#6BAED6", "#4292C6",
                                        "#2171B5", "#084594"))(100),
             main = "Subgenome Dominance by Tissue\n(proportion of group expression)",
             fontsize_row = 5,
             border_color = NA,
             na_col = "grey90")
    dev.off()
    message("#     → dominance_heatmap.pdf")
}, error = function(e) message(paste0("#     ERROR: ", e$message)))

# ---- 4. Collapsed PCA (Level 2) ---------------------------------------------
message("#   4/5 Collapsed PCA...")

# PCA on log2(TPM+1) of collapsed groups
log_collapsed_tpm <- log2(collapsed_tpm + 1)

# Only use groups with some expression
expressed_mask <- rowSums(collapsed_tpm) > 0
log_collapsed_expr <- log_collapsed_tpm[expressed_mask, ]

if (nrow(log_collapsed_expr) >= 3) {
    pca_collapsed <- prcomp(t(log_collapsed_expr), scale. = TRUE)
    pca_df <- data.frame(
        PC1 = pca_collapsed$x[, 1],
        PC2 = pca_collapsed$x[, 2],
        sample = rownames(pca_collapsed$x),
        stringsAsFactors = FALSE
    )
    pca_df$tissue <- sample_tissue[pca_df$sample]

    var_explained <- round(summary(pca_collapsed)$importance[2, 1:2] * 100, 1)

    tryCatch({
        p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = tissue)) +
            geom_point(size = 3.5) +
            labs(title = "PCA — Homeolog Group Expression (Collapsed Level)",
                 subtitle = "Based on log2(TPM+1) of summed homeolog groups",
                 x = paste0("PC1 (", var_explained[1], "% variance)"),
                 y = paste0("PC2 (", var_explained[2], "% variance)")) +
            theme_minimal() +
            theme(legend.position = "right")

        ggsave(file.path(RESULTS_DIR, "homeolog_pca.pdf"), p_pca,
               width = 8, height = 6)
        message("#     → homeolog_pca.pdf")
    }, error = function(e) message(paste0("#     ERROR: ", e$message)))
}

# ---- 5. Global Subgenome Dominance Summary -----------------------------------
message("#   5/5 Global subgenome summary...")

# Pie/bar chart of how many groups each subgenome dominates
global_dom <- dominant_subgenome %>%
    group_by(dominant_subgenome) %>%
    summarise(n_groups = n(), .groups = "drop") %>%
    mutate(fraction = n_groups / sum(n_groups))

tryCatch({
    p_global <- ggplot(global_dom,
                       aes(x = dominant_subgenome, y = n_groups,
                           fill = dominant_subgenome)) +
        geom_bar(stat = "identity", width = 0.6) +
        geom_text(aes(label = paste0(n_groups, " groups\n(",
                                     round(fraction * 100, 0), "%)")),
                  vjust = -0.3, size = 3.5) +
        scale_fill_manual(values = SUBGENOME_COLORS) +
        labs(title = "Subgenome Dominance Summary",
             subtitle = "Number of homeolog groups where each subgenome has highest expression",
             x = "Dominant Subgenome", y = "Number of Homeolog Groups") +
        theme_minimal() +
        theme(legend.position = "none") +
        ylim(0, max(global_dom$n_groups) * 1.2)

    ggsave(file.path(RESULTS_DIR, "subgenome_summary.pdf"), p_global,
           width = 6, height = 5)
    message("#     → subgenome_summary.pdf")
}, error = function(e) message(paste0("#     ERROR: ", e$message)))

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
message("\n# ── Summary Statistics ──")

# Genes expressed (TPM > 1 in at least one tissue)
tpm_aqp_mean <- data.frame(row.names = rownames(tpm_aqp))
for (tissue in unique(design$tissue)) {
    samps <- design$sample[design$tissue == tissue]
    samps <- samps[samps %in% colnames(tpm_aqp)]
    if (length(samps) > 0) {
        tpm_aqp_mean[[tissue]] <- rowMeans(tpm_aqp[, samps, drop = FALSE])
    }
}
n_expressed <- sum(apply(tpm_aqp_mean, 1, max) > 1)
message(paste0("#   Individual genes expressed (TPM>1): ", n_expressed, "/",
               length(aqp_in_counts)))

# Collapsed groups expressed
tpm_collapsed_mean <- data.frame(row.names = rownames(collapsed_tpm))
for (tissue in unique(design$tissue)) {
    samps <- design$sample[design$tissue == tissue]
    samps <- samps[samps %in% colnames(collapsed_tpm)]
    if (length(samps) > 0) {
        tpm_collapsed_mean[[tissue]] <- rowMeans(
            collapsed_tpm[, samps, drop = FALSE])
    }
}
n_groups_expressed <- sum(apply(tpm_collapsed_mean, 1, max) > 1)
message(paste0("#   Collapsed groups expressed (TPM>1): ", n_groups_expressed, "/",
               nrow(collapsed_tpm)))

# Dominance patterns
message(paste0("#   Subgenome dominance distribution:"))
for (sg in c("A", "B", "C", "D")) {
    n_dom <- sum(dominant_subgenome$dominant_subgenome == sg, na.rm = TRUE)
    message(paste0("#     Subgenome ", sg, ": dominant in ", n_dom, " groups"))
}

# Mean dominance proportion of the leading subgenome
mean_dom_prop <- mean(dominant_subgenome$dominance_proportion, na.rm = TRUE)
message(paste0("#   Mean dominance proportion: ",
               round(mean_dom_prop * 100, 1), "%"))
message(paste0("#   (25% = equal, >50% = strong dominance)"))

# Save combined summary
summary_stats <- data.frame(
    metric = c("total_aquaporin_genes", "genes_in_counts",
               "homeolog_groups", "complete_quartets",
               "genes_expressed_tpm1", "groups_expressed_tpm1",
               "mean_dominant_proportion",
               paste0("dominant_", c("A", "B", "C", "D"))),
    value = c(nrow(hg), length(aqp_in_counts),
              length(unique(hg$homeolog_group)),
              sum(group_annotation$completeness == "complete"),
              n_expressed, n_groups_expressed,
              round(mean_dom_prop, 4),
              sapply(c("A","B","C","D"), function(sg)
                  sum(dominant_subgenome$dominant_subgenome == sg, na.rm = TRUE)))
)
write.csv(summary_stats, file.path(RESULTS_DIR, "summary_statistics.csv"),
          row.names = FALSE)

message("\n# =================================================================")
message("# HOMEOLOG EXPRESSION ANALYSIS COMPLETE")
message(paste0("# Output directory: ", RESULTS_DIR))
message("# =================================================================")
