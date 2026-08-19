############################################################
# 02. Gene Annotation and Differentially Expressed Genes
############################################################


# ----------------------------------------------------------
# 02.1 Load DESeq2 results
# ----------------------------------------------------------

# Load the DESeq2 object generated in the previous step.

dds <- readRDS(
  "results/dds.rds"
)


# ----------------------------------------------------------
# 02.2 Define differential expression contrast
# ----------------------------------------------------------

# Compare the disease/case group against the control group.
# The reference level is Control.

res <- results(
  dds,
  contrast = c(
    "condition",
    "Case",
    "Control"
  )
)


# ----------------------------------------------------------
# 02.3 Convert results to data frame
# ----------------------------------------------------------

res_df <- as.data.frame(res)

res_df$Gene_ID <- rownames(res_df)

res_df <- res_df %>%
  dplyr::select(
    Gene_ID,
    everything()
  )


# ----------------------------------------------------------
# 02.4 Remove genes without statistical results
# ----------------------------------------------------------

res_df <- res_df %>%
  filter(
    !is.na(padj)
  )


# ----------------------------------------------------------
# 02.5 Define differentially expressed genes
# ----------------------------------------------------------

# Differentially expressed genes are defined using:
# |log2FoldChange| >= 1
# adjusted P-value < 0.05

DEGs <- res_df %>%
  filter(
    abs(log2FoldChange) >= 1,
    padj < 0.05
  )


# ----------------------------------------------------------
# 02.6 Classify DEGs by expression direction
# ----------------------------------------------------------

DEGs <- DEGs %>%
  mutate(
    Regulation = case_when(
      log2FoldChange >= 1 ~ "Upregulated",
      log2FoldChange <= -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )


# ----------------------------------------------------------
# 02.7 Extract upregulated genes
# ----------------------------------------------------------

upregulated_genes <- DEGs %>%
  filter(
    Regulation == "Upregulated"
  )


# ----------------------------------------------------------
# 02.8 Extract downregulated genes
# ----------------------------------------------------------

downregulated_genes <- DEGs %>%
  filter(
    Regulation == "Downregulated"
  )


# ----------------------------------------------------------
# 02.9 Summarize DEG results
# ----------------------------------------------------------

deg_summary <- data.frame(
  Category = c(
    "Total DEGs",
    "Upregulated",
    "Downregulated"
  ),
  Number = c(
    nrow(DEGs),
    nrow(upregulated_genes),
    nrow(downregulated_genes)
  )
)

print(deg_summary)


# ----------------------------------------------------------
# 02.10 Export complete DEG table
# ----------------------------------------------------------

dir.create(
  "results/DEGs",
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  DEGs,
  file = "results/DEGs/all_DEGs.csv",
  row.names = FALSE
)


# ----------------------------------------------------------
# 02.11 Export upregulated genes
# ----------------------------------------------------------

write.csv(
  upregulated_genes,
  file = "results/DEGs/upregulated_genes.csv",
  row.names = FALSE
)


# ----------------------------------------------------------
# 02.12 Export downregulated genes
# ----------------------------------------------------------

write.csv(
  downregulated_genes,
  file = "results/DEGs/downregulated_genes.csv",
  row.names = FALSE
)


# ----------------------------------------------------------
# 02.13 Export DEG summary
# ----------------------------------------------------------

write.csv(
  deg_summary,
  file = "results/DEGs/DEG_summary.csv",
  row.names = FALSE
)


# ----------------------------------------------------------
# 02.14 Prepare gene list for downstream analyses
# ----------------------------------------------------------

DEG_gene_list <- DEGs$Gene_ID

up_gene_list <- upregulated_genes$Gene_ID

down_gene_list <- downregulated_genes$Gene_ID


# ----------------------------------------------------------
# 02.15 Save DEG objects
# ----------------------------------------------------------

saveRDS(
  DEGs,
  file = "results/DEGs/DEGs.rds"
)

saveRDS(
  up_gene_list,
  file = "results/DEGs/upregulated_gene_list.rds"
)

saveRDS(
  down_gene_list,
  file = "results/DEGs/downregulated_gene_list.rds"
)
