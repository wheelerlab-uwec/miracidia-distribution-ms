library(tidyverse)
library(broom)
library(ggrepel)
library(DESeq2)
library(ggtext)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)

# keep arg for pipeline
input_file <- args[1]

# use below input file if starting with counts
# input_file <- here('pipeline', 'counts', 'star', 'counts.tsv')

counts_df <- read_tsv(input_file, comment = "#") |>
  mutate(across(where(is.numeric), as.integer))
counts_df

counts_summary <- counts_df |>
  select(Geneid, contains('bam')) |>
  rename_with(~ str_remove(., "dedup/star/.*.bam:"), everything()) |>
  rowwise() |>
  mutate(total_counts = sum(c_across(where(is.numeric)), na.rm = T)) |>
  filter(total_counts >= 10)
counts_summary

sample_summary <- counts_df |>
  select(Geneid, contains('bam')) |>
  rename_with(~ str_remove(., "dedup/star/.*.bam:"), everything()) |>
  pivot_longer(-Geneid, names_to = 'sample', values_to = 'count') |>
  filter(count > 0) |>
  group_by(Geneid) |>
  tally() |>
  filter(n <= 3)
sample_summary

genes_to_remove = sample_summary$Geneid

counts_filt <- counts_summary |>
  filter(!Geneid %in% genes_to_remove) |>
  arrange(Geneid) |>
  select(-total_counts)
counts_filt

counts_m <- counts_filt |>
  select(-Geneid) |>
  as.matrix()
rownames(counts_m) <- counts_filt$Geneid
counts_m

metadata <- data.frame(sample_id = colnames(counts_m)) |>
  mutate(tissue = str_sub(sample_id, 1, 3), rep = str_sub(sample_id, 6, 6))
rownames(metadata) <- metadata$sample_id
metadata <- select(metadata, -sample_id)

dds <- DESeqDataSetFromMatrix(countData = counts_m, colData = metadata, design = ~tissue)
dds <- DESeq(dds)

res <- results(dds)


###### PCA

vsd <- vst(dds)

pca_fit <- t(assay(vsd)) |>
  prcomp(scale = TRUE)
# pca_fit

(pca_plot <- pca_fit |>
  augment(t(assay(vsd))) |>
  dplyr::rename(sample = .rownames) |>
  separate(sample, into = c('tissue', 'rep'), sep = '-') |>
  mutate(
    tissue = case_when(
      tissue == 'Int' ~ 'Intestine',
      tissue == 'Liv' ~ 'Liver'
    ),
    rep = str_remove(rep, "_.*")
  ) |>
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, color = tissue)) +
  geom_point(aes(shape = as.factor(rep)), size = 4, show.legend = FALSE) +
  stat_ellipse(type = "norm", level = 0.68) + # 1 SD
  scale_color_manual(values = c('indianred', 'steelblue')) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  labs(x = "PC1", y = "PC2", color = "Tissue source", shape = "Replicate") +
  cowplot::theme_half_open() +
  theme(
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_markdown(size = 8),
    legend.position = 'bottom'
  ) +
  NULL)

ggsave('pipeline/plots/pca.pdf', pca_plot, width = 6, height = 6)
ggsave('pipeline/plots/pca.png', pca_plot, width = 6, height = 6, bg = 'white')

###### Diff Exp

volcano_data <- as_tibble(res, rownames = 'gene_id')

degs <- volcano_data |>
  filter(abs(log2FoldChange) > 1, padj < 0.05)

(volcano_plot <- volcano_data |>
  mutate(
    color = case_when(
      log2FoldChange > 1 & padj < 0.05 | log2FoldChange < -1 & padj < 0.05 ~ 'deg',
      TRUE ~ 'not deg'
    )
  ) |>
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = color), alpha = 0.5, show.legend = FALSE) +
  # geom_label_repel(data = degs,
  #                  aes(x = log2FoldChange, y = -log10(padj), label = gene_id),
  #                  max.overlaps = 50, size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_vline(xintercept = -1, linetype = 'dashed') +
  labs(x = "log<sub>2</sub>(Fold change)", y = '-log<sub>10</sub>(Adjusted p-value)') +
  scale_color_manual(values = c('#EBB6B3', '#334139')) +
  cowplot::theme_half_open() +
  theme(
    axis.title = element_markdown(size = 9),
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    legend.position = 'bottom'
  ) +
  NULL)

ggsave('pipeline/plots/volcano.png', volcano_plot, bg = 'white', width = 6, height = 6)

write_csv(volcano_data, here('pipeline/deseq_results/deseq_results.csv'))

norm_counts <- counts(dds, normalized = TRUE) |>
  as_tibble(rownames = 'gene_id')
write_csv(norm_counts, here('pipeline/deseq_results/deseq_norm_counts.csv'))

fig1 <- plot_grid(
  pca_plot,
  volcano_plot,
  align = 'h',
  axis = 'tb',
  labels = c('A', 'B')
)

save_plot(here('Fig1/Fig1.pdf'), fig1, base_width = 6, base_height = 3)

save_plot(here('Fig1/Fig1.png'), fig1, base_width = 6, base_height = 3, bg = 'white')
###### comparison to eggs

egg_degs <- read_tsv('Fig1/egg_degs.tsv') |>
  select(gene_id = Gene, padj, log2FoldChange) |>
  filter(abs(log2FoldChange) > 1, padj < 0.05)


mira_degs <- degs |>
  mutate(
    shared = case_when(
      gene_id %in% egg_degs$gene_id ~ TRUE,
      TRUE ~ FALSE
    )
  ) |>
  write_csv("Fig1/Table1.csv")
