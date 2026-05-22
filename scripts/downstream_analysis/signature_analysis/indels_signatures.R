#!/usr/bin/env Rscript

# ==============================================================================
# 1. ESSENTIAL PACKAGES
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(yaml)
  library(rstudioapi)
  library(scales)
  
  # Genomics & Mutational Signatures
  library(MutationalPatterns)
  library(GenomicRanges)
  library(BSgenome)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

# ==============================================================================
# 2. DIRECTORY CONFIGURATION & HELPERS
# ==============================================================================
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
}

config     <- yaml::read_yaml(file.path(script_dir, "config-indel.yaml"))
vcf_dir    <- file.path(script_dir, config$vcf_dir)
group_file <- file.path(vcf_dir, config$indel_file)

cat("[INFO] VCF dir:", vcf_dir, "\n")
cat("[INFO] Group file:", group_file, "\n")

outdir     <- file.path(script_dir, "results")
plots_dir  <- file.path(outdir, "plots")
tables_dir <- file.path(outdir, "tables")
rds_dir    <- file.path(outdir, "rds")

invisible(lapply(c(plots_dir, tables_dir, rds_dir), dir.create, showWarnings = FALSE, recursive = TRUE))

save_table <- function(df, name) fwrite(df, file.path(tables_dir, paste0(name, ".tsv")), sep = "\t")
save_rds   <- function(obj, name) saveRDS(obj, file.path(rds_dir, paste0(name, ".rds")))
save_plot  <- function(plot, name, width = 10, height = 10) {
  ggsave(filename = file.path(plots_dir, paste0(name, ".png")), plot = plot, width = width, height = height, dpi = 300)
}

# ==============================================================================
# 3. GENOMIC REFERENCE & DATA LOADING
# ==============================================================================
seqinfo_hg38 <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

group_info <- fread(group_file, header = FALSE, col.names = c("groups", "sample"))

samples_per_group <- group_info %>%
  group_by(groups) %>%
  summarise(count = n_distinct(sample), .groups = "drop")

sample_names <- group_info$sample

# --- NAME CLEANING ---
sample_id <- gsub("indels.hg38_multianno.txt", "", sample_names)
sample_id <- gsub("\\.txt$", "", sample_id)
sample_id <- gsub("\\.tsv$", "", sample_id)

groups <- group_info$groups

myfiles <- lapply(file.path(vcf_dir, sample_names), function(x) {
  fread(x, header = TRUE, sep = "\t", colClasses = c("Ref" = "character", "Alt" = "character"))
})

# ==============================================================================
# 4. GRANGES LIST CREATION
# ==============================================================================
# --- PER SAMPLE (Maintained for logging and QC) ---
sample_grList <- GRangesList(lapply(seq_along(sample_names), function(i) {
  df <- as.data.frame(myfiles[[i]])[, 1:5]
  gr <- GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), ref = df[[4]], alt = df[[5]], seqinfo = seqinfo_hg38)
  mcols(gr)$group    <- groups[i]
  mcols(gr)$sampleID <- sample_id[i]
  return(gr)
}))
names(sample_grList) <- sample_id
save_rds(sample_grList, "indel_sample_grList")

# --- PER GROUP (Combined list used for plotting) ---
unique_groups <- unique(groups)
groups_grList <- GRangesList(lapply(unique_groups, function(g) {
  idx <- which(groups == g)
  df  <- bind_rows(lapply(idx, function(i) {
    tmp <- as.data.frame(myfiles[[i]])[, 1:5]
    tmp$sampleID <- sample_id[i]
    tmp$group <- groups[i]
    tmp
  }))
  gr  <- GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), ref = df[[4]], alt = df[[5]], seqinfo = seqinfo_hg38)
  mcols(gr)$sampleID <- df$sampleID
  mcols(gr)$group    <- df$group
  return(gr)
}))
names(groups_grList) <- as.vector(unique_groups)
save_rds(groups_grList, "indel_groups_grList")

# ==============================================================================
# 5. QC: RECORD COUNT CHECK
# ==============================================================================
sample_summary <- tibble(
  sampleID = names(sample_grList),
  group = sapply(sample_grList, function(gr) unique(mcols(gr)$group)),
  records = lengths(sample_grList)
)

group_summary <- tibble(group = names(groups_grList), records = lengths(groups_grList))

qc_check <- group_summary %>%
  left_join(sample_summary %>% group_by(group) %>% summarise(calc_records = sum(records)), by = "group") %>%
  mutate(match = records == calc_records, status = ifelse(match, "OK", "Error"))

save_table(qc_check, "indel_qc_check")
cat("\n[INFO] Final QC Match Status:", qc_check$status, "\n")

# ==============================================================================
# 6. MUTATIONAL PATTERNS: INDEL CONTEXTS
# ==============================================================================
cat("\n[INFO] Generating Indel Contexts...\n")

# Using groups_grList so counts matrix columns represent GROUPS
indel_grl    <- get_indel_context(get_mut_type(groups_grList, type = "indel"), "BSgenome.Hsapiens.UCSC.hg38")
indel_counts <- count_indel_contexts(indel_grl)

# Built-in Plots
plot_indel <- plot_indel_contexts(indel_counts, condensed = FALSE, same_y = TRUE) +
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 90, hjust = 0.5), axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 3, face = "bold"), axis.title.y = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
        strip.text.x = element_text(size = 20, face = "plain", color = "black"),
        strip.text.y = element_text(size = 20, face = "plain", color = "black"))
save_plot(plot_indel, "indel_condensed", width = 25, height = 7)

plot_main_contexts <- plot_main_indel_contexts(indel_counts, same_y = TRUE) +
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 15, hjust = 0.5), axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 3, face = "bold"), axis.title.y = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
        strip.text.x = element_text(size = 20, face = "bold", color = "black"),
        strip.text.y = element_text(size = 20, face = "bold", color = "black"))
save_plot(plot_main_contexts, "indel_main_context", height = 10, width = 15)

# ==============================================================================
# 7. COSMIC SIGNATURE REFITTING (Run Once)
# ==============================================================================
cat("\n[INFO] Fitting COSMIC Indel Signatures...\n")

signatures_indel <- get_known_signatures(muttype = "indel")

# Broad fit to identify candidate signatures
fit_res_indel <- fit_to_signatures(indel_counts, signatures_indel)
save_rds(fit_res_indel, "fit_res_indel")

selected_sigs_indel <- which(rowSums(fit_res_indel$contribution) > 0)
signatures_indel_reduced <- signatures_indel[, selected_sigs_indel]

# Strict refit for higher precision
strict_refit_indel <- fit_to_signatures_strict(
  mut_mat = indel_counts,
  signatures = signatures_indel_reduced,
  max_delta = 0.03,
  method = "best_subset"
)
fit_res_indel_strict <- strict_refit_indel$fit_res

# ==============================================================================
# 8. CUSTOM PLOTTING (Transparency, Borders, & Labels)
# ==============================================================================
cat("\n[INFO] Generating Custom Publication Plots...\n")

# Prepare Data for ggplot
contrib_long <- as.data.frame(fit_res_indel_strict$contribution) %>%
  tibble::rownames_to_column("Signature") %>%
  pivot_longer(cols = -Signature, names_to = "Group", values_to = "Absolute") %>%
  filter(Absolute > 0) %>% # Clear legend: only signatures with > 0 contribution
  group_by(Group) %>%
  mutate(Total = sum(Absolute),
         Relative = Absolute / Total,
         Label = ifelse(Relative > 0.04, paste0(round(Relative * 100, 1), "%"), "")) %>%
  ungroup()

# Force signature order for consistent colors
ordered_sigs <- unique(contrib_long$Signature)
contrib_long$Signature <- factor(contrib_long$Signature, levels = ordered_sigs)

# Custom Indel Palette
my_colors_indel <- c("#F39C12FF", "#2471a3", "#d98880", "#7E6148FF", "#B09C85FF", 
                     "#8491B4FF", "#91D1C2FF", "#CD534CFF")
final_palette <- setNames(my_colors_indel[seq_along(ordered_sigs)], ordered_sigs)

theme_indel_clean <- theme_classic() + theme(
  legend.position = "right",
  legend.title = element_text(size = 18, face = "bold"),
  legend.text = element_text(size = 16),
  axis.text.x = element_text(size = 20, color = "black", face = "bold"),
  axis.text.y = element_text(size = 16, color = "black"),
  axis.title = element_text(size = 18, face = "bold"),
  panel.grid.major.y = element_line(color = "gray90", linetype = "dotted"),
  plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
)

# --- RELATIVE Plot (%) ---
p_rel_clean <- ggplot(contrib_long, aes(x = Group, y = Relative, fill = Signature)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.6, alpha = 0.7) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 6, color = "black", fontface = "bold") +
  scale_fill_manual(values = final_palette, name = "ID COSMIC") +
  scale_y_continuous(labels = percent) +
  labs(title = "", x = "", y = "Relative Contribution") + theme_indel_clean

# --- ABSOLUTE Plot (Counts) ---
p_abs_clean <- ggplot(contrib_long, aes(x = Group, y = Absolute, fill = Signature)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.6, alpha = 0.7) +
  scale_fill_manual(values = final_palette, name = "ID COSMIC") +
  labs(title = "", x = "", y = "Absolute Contribution") + theme_indel_clean

# --- Combine with Patchwork ---
combined_indel <- (p_rel_clean | p_abs_clean) + 
  plot_annotation(
    title = "Indel Mutational Signatures (COSMIC v3.5)",
    subtitle = NULL,
    theme = theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 14, hjust = 0.5))
  ) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "right", legend.key.size = unit(1.2, "cm"))

print(combined_indel)
save_plot(combined_indel, "indel_contribution_FINAL_CLEAN", width = 14, height = 7)

cat("\n[INFO] Script finished successfully.\n")
