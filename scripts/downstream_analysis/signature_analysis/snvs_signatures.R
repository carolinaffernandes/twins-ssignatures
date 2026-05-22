#!/usr/bin/env Rscript

# ==============================================================================
# 1. PACKAGES
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
  library(readxl)
  library(yaml)
  library(rstudioapi)
  
  # Visualization
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(gridExtra)
  library(RColorBrewer)
  library(ggrepel)
  
  # Genomics & Mutational Signatures
  library(MutationalPatterns)
  library(SomaticSignatures)
  library(GenomicRanges)
  library(VariantAnnotation)
  library(BSgenome)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(NMF)
})

# ==============================================================================
# 2. DIRECTORY CONFIGURATION & HELPERS
# ==============================================================================
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
}

config     <- yaml::read_yaml(file.path(script_dir, "config.yaml"))
vcf_dir    <- file.path(script_dir, config$vcf_dir)
group_file <- file.path(vcf_dir, config$group_file)

cat("[INFO] Script dir:", script_dir, "\n[INFO] VCF dir:", vcf_dir, "\n")

outdir     <- file.path(script_dir, "results")
plots_dir  <- file.path(outdir, "plots")
tables_dir <- file.path(outdir, "tables")
rds_dir    <- file.path(outdir, "rds")

invisible(lapply(c(plots_dir, tables_dir, rds_dir), dir.create, showWarnings = FALSE, recursive = TRUE))

save_table <- function(df, name) fwrite(df, file.path(tables_dir, paste0(name, ".tsv")), sep = "\t")
save_rds   <- function(obj, name) saveRDS(obj, file.path(rds_dir, paste0(name, ".rds")))
save_plot  <- function(plot, name, width = 10, height = 10) {
  ggsave(filename = file.path(plots_dir, paste0(name, ".png")), plot = plot, width = width, height = height)
}

# ==============================================================================
# 3. DATA LOADING & PREPARATION
# ==============================================================================
seqinfo_hg38 <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

# Group Metadata
group_info <- fread(group_file, header = FALSE, col.names = c("groups", "sample"))
sample_id  <- sub(".snps.hg38_multianno.txt", "", group_info$sample)
groups     <- group_info$groups

cat("[INFO] Samples found:", length(sample_id), "\n")

# Read variants
myfiles <- lapply(file.path(vcf_dir, group_info$sample), function(x) {
  fread(x, header = TRUE, sep = "\t", colClasses = c("Ref" = "character", "Alt" = "character"))
})

# Create GRanges (By Sample)
sample_grList <- GRangesList(lapply(seq_along(sample_id), function(i) {
  df <- as.data.frame(myfiles[[i]])[, 1:5]
  gr <- GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), ref = df[[4]], alt = df[[5]], seqinfo = seqinfo_hg38)
  mcols(gr)$group <- groups[i]
  mcols(gr)$sampleID <- sample_id[i]
  return(gr)
}))
names(sample_grList) <- sample_id
save_rds(sample_grList, "sample_grList")

# Create GRanges (By Group)
unique_groups <- unique(groups)
groups_grList <- GRangesList(lapply(unique_groups, function(g) {
  idx <- which(groups == g)
  df  <- bind_rows(lapply(idx, function(i) as.data.frame(myfiles[[i]])[, 1:5]))
  gr  <- GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), ref = df[[4]], alt = df[[5]], seqinfo = seqinfo_hg38)
  mcols(gr)$group <- g
  return(gr)
}))
names(groups_grList) <- unique_groups
save_rds(groups_grList, "groups_grList")

# ==============================================================================
# 4. MUTATIONAL PATTERNS: SPECTRUM & MATRIX
# ==============================================================================
cat("\n[INFO] Generating Mutational Matrices and Spectra...\n")

# By Group
type_occ_groups <- mut_type_occurrences(groups_grList, "BSgenome.Hsapiens.UCSC.hg38")
save_plot(plot_spectrum(type_occ_groups, CT = TRUE, by = unique_groups, legend = TRUE), "snv_spectrum_by_group", 12, 8)

mut_mat <- mut_matrix(vcf_list = groups_grList, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
colnames(mut_mat) <- unique_groups 
save_plot(plot_96_profile(mut_mat, ymax = 0.08), "snv_profile_96_ct_GROUPS", 15, 10)

# ==============================================================================
# 5. COSMIC SIGNATURE REFIT (V3.5)
# ==============================================================================
cat("\n[INFO] Fitting COSMIC Signatures...\n")

path_cosmic <- "/home/carolis/COSMIC_catalogue-signatures_SBS96_v3.5/COSMIC_v3.5_SBS_GRCh38.txt"
cosmic_v35  <- read.table(path_cosmic, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
cosmic_v35  <- as.matrix(cosmic_v35)[rownames(mut_mat), ]

# Fit
fit_strict  <- fit_to_signatures_strict(mut_mat, cosmic_v35, max_delta = 0.048)$fit_res
colnames(fit_strict$contribution) <- unique_groups
selected_sig <- which(rowSums(fit_strict$contribution) > 0)

# Save Raw Data
write.table(fit_strict$contribution[selected_sig, , drop=FALSE],
            file.path(tables_dir, "snv_mutationalSignaturesABSOLUTE_groups.txt"),
            sep = "\t", row.names = TRUE, quote = FALSE)

# --- THEME & DATA PREP FOR COSMIC PLOTS ---
my_colors <- c("#D44E00", "#91D1C2FF", "gray80", "gray30", "#4DBBD5FF", "#00A087FF",
               "#868686FF", "#91D1C2FF", "#8491B4FF", "#7E6148FF", "#B09C85FF", "#0073C2FF", "#868686FF")

theme_clean <- theme_classic(base_size = 16) + 
  theme(
    legend.position = "right", legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14), axis.text.x = element_text(size = 20, face = "bold", color = "black"), 
    axis.text.y = element_text(size = 18, color = "black"), axis.title = element_text(size = 18, face = "bold"),
    axis.line = element_line(linewidth = 0.8), panel.grid.major.y = element_line(color = "gray95", linetype = "dotted"), 
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5), plot.margin = margin(10, 10, 10, 10)
  )

contrib_data <- as.data.frame(fit_strict$contribution[selected_sig, , drop=FALSE]) %>% rownames_to_column("Signature")
ordered_sigs <- contrib_data$Signature
final_palette <- setNames(my_colors[seq_along(ordered_sigs)], ordered_sigs)

contrib_long <- contrib_data %>%
  pivot_longer(cols = -Signature, names_to = "Sample", values_to = "Absolute") %>%
  mutate(Signature = factor(Signature, levels = ordered_sigs)) %>%
  group_by(Sample) %>%
  mutate(Relative = Absolute / sum(Absolute),
         Label = ifelse(Relative > 0.03, paste0(round(Relative * 100, 1), "%"), "")) %>%
  ungroup()

# --- PLOT 1: SOLID BARS (Percentages) ---
p_rel_solid <- ggplot(contrib_long, aes(x = Sample, y = Relative, fill = Signature)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 6, fontface = "bold") +
  scale_fill_manual(values = final_palette) + scale_y_continuous(labels = percent, expand = expansion(c(0, 0.05))) +
  labs(x = NULL, y = "Relative Contribution") + theme_clean

p_abs_solid <- ggplot(contrib_long, aes(x = Sample, y = Absolute, fill = Signature)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = final_palette) + scale_y_continuous(expand = expansion(c(0, 0.1))) +
  labs(x = NULL, y = "Absolute Contribution") + theme_clean

combo_solid <- (p_rel_solid | p_abs_solid) + plot_annotation(title = "Mutational Signatures Analysis (COSMIC v3.5)")
ggsave(file.path(plots_dir, "SNV_Signatures_Final_Percentages.png"), combo_solid, width = 18, height = 9, dpi = 300)

# --- PLOT 2: TRANSPARENT BARS ---
p_rel_trans <- ggplot(contrib_long, aes(x = Sample, y = Relative, fill = Signature)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.6, alpha = 0.7) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 6, fontface = "bold") +
  scale_fill_manual(values = final_palette) + scale_y_continuous(labels = percent, expand = expansion(c(0, 0.05))) +
  labs(x = NULL, y = "Relative Contribution") + theme_clean

p_abs_trans <- ggplot(contrib_long, aes(x = Sample, y = Absolute, fill = Signature)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.6, alpha = 0.7) +
  scale_fill_manual(values = final_palette) + scale_y_continuous(expand = expansion(c(0, 0.1))) +
  labs(x = NULL, y = "Absolute Contribution") + theme_clean

combo_trans <- (p_rel_trans | p_abs_trans) + 
  plot_annotation(title = "SBS Mutational Signature (COSMIC v3.5)") +
  theme(legend.key.size = unit(1.2, "cm"))
ggsave(file.path(plots_dir, "SNV_Signatures_Final_TRANSPARENCY.png"), combo_trans, width = 14, height = 7, dpi = 300)

# ==============================================================================
# 6. SOMATIC SIGNATURES (NMF & OBSERVED SPECTRUM)
# ==============================================================================
cat("\n[INFO] Running Somatic Signatures (NMF)...\n")

# Prepare VRanges from existing 'myfiles'
merged_df <- bind_rows(lapply(seq_along(sample_id), function(i) {
  df <- as.data.frame(myfiles[[i]])[, 1:5]
  colnames(df) <- c("chr", "startvar", "endvar", "ref", "var")
  df$sample_id <- sample_id[i]
  df$group <- groups[i]
  return(df)
}))

merged_vr <- VRanges(seqnames = merged_df$chr, ranges = IRanges(merged_df$startvar, merged_df$endvar),
                     ref = merged_df$ref, alt = merged_df$var, study = merged_df$group, sampleNames = merged_df$sample_id)

merged_motifs <- mutationContext(merged_vr, BSgenome.Hsapiens.UCSC.hg38)
merged_mm     <- motifMatrix(merged_motifs, group = "study", normalize = TRUE)
merged_mm     <- merged_mm[rowSums(merged_mm == 0) != ncol(merged_mm), ]

sigs_nmf <- identifySignatures(merged_mm, 2, nmfDecomposition)

# Custom Observed Spectrum Plot
p_nmf <- plotObservedSpectrum(sigs_nmf) + 
  facet_grid(factor(sample, levels = c('M3', 'M6')) ~ alteration) + scale_y_continuous(limits = c(0, 0.04))

p_nmf$layers[[1]]$aes_params$colour <- "black"
p_nmf$layers[[1]]$aes_params$linewidth <- 0.6 
p_nmf$layers[[1]]$aes_params$alpha <- 0.7

p_nmf <- p_nmf + labs(x = NULL, y = "Contribution") +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        axis.ticks = element_line(color = 'black', linewidth = 0.2),
        axis.title.y = element_text(size = 12, color = 'black', face = "bold"),
        axis.text.y = element_text(size = 10, color = 'black'), axis.text.x = element_text(size = 8, color = 'black'),
        strip.text.y = element_text(size = 12, face = 'bold'), strip.background.y = element_rect(color = "black", fill = "gray"),
        strip.text.x = element_text(size = 12, face = 'plain')) + 
  scale_fill_manual(values = c("M3" = "#0072B2", "M6" = "#D55E00"))

ggsave(file.path(plots_dir, "SomaticSig_ObservedSpectrum_SNV_Customized.png"), p_nmf, dpi = 600, width = 30, height = 12, units = 'cm')

# ==============================================================================
# 7. COSINE SIMILARITY
# ==============================================================================
cat("\n[INFO] Calculating Cosine Similarity...\n")

df_corr <- as.data.frame(motifMatrix(merged_motifs, group = "study", normalize = FALSE)) %>%
  mutate(M3_prop = M3 / sum(M3), M6_prop = M6 / sum(M6), 
         diff_abs = abs(M3_prop - M6_prop), Significant = diff_abs > 0.015, Motif = rownames(.))

cos_sim <- sum(df_corr$M3_prop * df_corr$M6_prop) / (sqrt(sum(df_corr$M3_prop^2)) * sqrt(sum(df_corr$M6_prop^2)))

plot_cosine <- ggplot(df_corr, aes(x = M3_prop, y = M6_prop)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray70") +
  geom_point(aes(color = Significant), alpha = 0.6, size = 2.5) +
  scale_color_manual(values = c("FALSE" = "#2c3e50", "TRUE" = "#e74c3c")) +
  geom_text_repel(data = subset(df_corr, Significant == TRUE), aes(label = Motif), size = 3.5, fontface = "bold", box.padding = 0.5) +
  labs(title = "SBS96 Mutational Profile Conservation", subtitle = sprintf("Cosine Similarity = %.4f", cos_sim),
       x = "Relative Frequency in M3", y = "Relative Frequency in M6") +
  theme_classic() + theme(legend.position = "none", panel.grid = element_blank())

ggsave(file.path(plots_dir, "M3_M6_Cosine_Correlation_Plot.png"), plot_cosine, width = 7, height = 7, dpi = 300)

cat("\n[INFO] Script finished successfully.\n")
