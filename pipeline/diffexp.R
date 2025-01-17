library(tidyverse)
library(broom)
library(ggrepel)
library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]

counts_df <- read_tsv(input_file,
                      comment = "#") |>
             mutate(across(where(is.numeric), as.integer))
counts_df

counts_summary <- counts_df |>
    select(Geneid, contains('bam')) |>
    rename_with(~str_remove(., "dedup/star/.*.bam:"), everything()) |>
    rowwise() |>
    mutate(total_counts = sum(c_across(where(is.numeric)), na.rm = T)) |>
    filter(total_counts >= 10)
counts_summary

sample_summary <- counts_df |>
    select(Geneid, contains('bam')) |>
    rename_with(~str_remove(., "dedup/star/.*.bam:"), everything()) |>
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
    mutate(tissue = str_sub(sample_id, 1, 3),
           rep = str_sub(sample_id, 6, 6))
rownames(metadata) <- metadata$sample_id
metadata <- select(metadata, -sample_id)

dds <- DESeqDataSetFromMatrix(countData = counts_m,
                              colData = metadata,
                              design = ~ tissue)
dds <- DESeq(dds)

res <- results(dds)



###### PCA

vsd <- vst(dds)

pca_fit <- t(assay(vsd)) |> 
  prcomp(scale = TRUE)
pca_fit

(pca_plot <- pca_fit |>
  augment(t(assay(vsd))) |>
  dplyr::rename(sample = .rownames) |>
  separate(sample, into = c('tissue', 'rep'), sep = '-') |>
    mutate(tissue = case_when(
    tissue == 'Int' ~ 'Intestine',
    tissue == 'Liv' ~ 'Liver'
  )) |>
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, color = tissue, shape = tissue)) + 
  geom_point(size = 4) +
  scale_color_manual(values = c('indianred', 'steelblue')) +
  labs(x = "PC1 (29% of variance)", y = "PC2 (23% of variance)", color = "Tissue source", shape = "Tissue source") +
  theme_minimal() +
  NULL
)

ggsave('plots/pca.pdf', pca_plot, width = 4, height = 4)
ggsave('plots/pca.png', pca_plot, width = 4, height = 4)

###### Diff Exp

volcano_data <- as_tibble(res, rownames = 'gene_id')
degs <- volcano_data |>
    filter(log2FoldChange > 2 | log2FoldChange < -2, 
           -log10(padj) > -log10(.05))

volcano_plot <- volcano_data |> 
    ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point() +
    geom_label_repel(data = degs,
                     aes(x = log2FoldChange, y = -log10(padj), label = gene_id),
                     max.overlaps = 50, size = 1) +
    geom_hline(yintercept = -log10(0.05)) +
    geom_vline(xintercept = 2) +
    geom_vline(xintercept = -2) +
    theme_minimal()

ggsave('plots/volcano.png', volcano_plot, bg = 'white', width = 6, height = 6)















# vsd_dists <- dist(t(assay(vsd)))

# vsd_dists_df <- as.matrix(vsd_dists) |>
#     as_tibble(rownames = 'sample')

# vsd_dist_plot <- vsd_dists_df |>
#     pivot_longer(-sample, names_to = 'comp', values_to = 'dist') |>
#     ggplot(aes(x = sample, y = comp, fill = dist)) +
#     geom_tile() +
#     scale_fill_viridis_c() +
#     coord_equal() +
#     NULL
# vsd_dist_plot
