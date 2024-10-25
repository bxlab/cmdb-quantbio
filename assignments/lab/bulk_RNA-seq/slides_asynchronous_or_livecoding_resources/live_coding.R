library(data.table)
library(tidyverse)
library(broom)
library(DESeq2)

# set your working directory to where your data and output will be stored
setwd("~/Dropbox/teaching/2024_qblab/differential_expression/")

# load the gene expression counts
counts_df <- read_delim("pasilla/pasilla_gene_counts.tsv")
# move the gene_id column to rownames, so that the contents of the
# tibble is entirely numeric
counts_df <- column_to_rownames(counts_df, var = "gene_id")

# look at first five rows
counts_df[1:5,]

# load the metadata
metadata_df <- read_delim("pasilla/pasilla_metadata.csv")
# move the sample IDs from the first column to rownames
metadata_df <- column_to_rownames(metadata_df, var = "SAMPLE_ID")

# check that the columns of the counts are identical and in the same order as the
# rows of the metadata
colnames(counts_df) == rownames(metadata_df)
table(colnames(counts_df) == rownames(metadata_df))

# create a DESeq object
dds <- DESeqDataSetFromMatrix(countData = counts_df,
                              colData = metadata_df,
                              design = ~ condition + type)

# apply VST normalization
vsd <- vst(dds)

# apply and plot principal components
plotPCA(vsd, intgroup = "condition")
plotPCA(vsd, intgroup = "type")
plotPCA(vsd, intgroup = c("condition", "type"))

# extract the normalized expression data and bind to metadata
vsd_df <- assay(vsd) %>%
  t() %>%
  as_tibble()
vsd_df <- bind_cols(metadata_df, vsd_df)

# examine the distribution of expression for a given gene
hist(vsd_df$FBgn0000008)

# apply multiple linear regression to a given gene
lm(data = vsd_df, formula = FBgn0000008 ~ condition + type) %>%
  summary() %>%
  tidy()

# use DESeq2 to perform differential expression analysis across all genes
dds <- DESeq(dds)

pasilla_res <- results(dds, name = "condition_untreated_vs_treated") %>%
  as_tibble(rownames = "GENE_ID")

pasilla_res <- pasilla_res %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

# if I wanted "untreated" to be the reference level
dds$condition <- relevel(dds$condition, ref = "untreated")

dds <- DESeq(dds)

pasilla_res <- results(dds, name = "condition_treated_vs_untreated") %>%
  as_tibble(rownames = "GENE_ID")

pasilla_res <- pasilla_res %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

# create a volcano plot

ggplot(data = pasilla_res, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = (abs(log2FoldChange) > 2 & pvalue < 1e-20))) +
  geom_text(data = pasilla_res %>% filter(abs(log2FoldChange) > 2 & pvalue < 1e-50),
            aes(x = log2FoldChange, y = -log10(pvalue) + 5, label = GENE_ID), size = 3,) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("darkgray", "coral")) +
  labs(y = expression(-log[10]("p-value")), x = expression(log[2]("fold change")))



