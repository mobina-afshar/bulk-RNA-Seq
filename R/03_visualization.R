############################################################
# 03. Differential Expression Visualization
############################################################


# ----------------------------------------------------------
# 03.1 Load required packages
# ----------------------------------------------------------

library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(grid)


# ----------------------------------------------------------
# 03.2 Prepare data for volcano plot
# ----------------------------------------------------------

volcano_data <- results_annotation %>%
  
  dplyr::mutate(
    
    padj_plot = ifelse(
      is.na(padj),
      1,
      padj
    ),
    
    padj_plot = ifelse(
      padj_plot == 0,
      .Machine$double.xmin,
      padj_plot
    ),
    
    negLog10Padj = -log10(padj_plot),
    
    Significance = dplyr::case_when(
      
      padj < 0.05 &
        log2FoldChange >= 1 ~ "Upregulated",
      
      padj < 0.05 &
        log2FoldChange <= -1 ~ "Downregulated",
      
      TRUE ~ "Not Significant"
      
    )
    
  )


# ----------------------------------------------------------
# 03.3 Select genes for labeling
# ----------------------------------------------------------

# Top 10 upregulated genes

top_up <- volcano_data %>%
  
  dplyr::filter(
    Significance == "Upregulated"
  ) %>%
  
  dplyr::arrange(padj) %>%
  
  dplyr::slice_head(n = 10)


# Top 10 downregulated genes

top_down <- volcano_data %>%
  
  dplyr::filter(
    Significance == "Downregulated"
  ) %>%
  
  dplyr::arrange(padj) %>%
  
  dplyr::slice_head(n = 10)


# Top 5 lncRNAs

top_lnc <- volcano_data %>%
  
  dplyr::filter(
    
    Gene_Biotype %in% c(
      
      "lncRNA",
      "antisense",
      "sense_intronic",
      "sense_overlapping",
      "processed_transcript",
      "bidirectional_promoter_lncRNA",
      "macro_lncRNA"
      
    ),
    
    padj < 0.05
    
  ) %>%
  
  dplyr::arrange(padj) %>%
  
  dplyr::slice_head(n = 5)


# ----------------------------------------------------------
# 03.4 Create final volcano plot labels
# ----------------------------------------------------------

label_data <- dplyr::bind_rows(
  
  top_up,
  top_down,
  top_lnc
  
) %>%
  
  dplyr::distinct(
    Gene_Symbol,
    .keep_all = TRUE
  )


# ----------------------------------------------------------
# 03.5 Generate volcano plot
# ----------------------------------------------------------

volcano_plot <- ggplot(
  
  volcano_data,
  
  aes(
    x = log2FoldChange,
    y = negLog10Padj
  )
  
) +
  
  geom_point(
    aes(
      color = Significance
    ),
    alpha = 0.75,
    size = 2
  ) +
  
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    colour = "grey60"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    colour = "grey60"
  ) +
  
  geom_text_repel(
    
    data = label_data,
    
    aes(
      label = Gene_Symbol
    ),
    
    size = 3.5,
    
    box.padding = 0.4,
    
    point.padding = 0.3,
    
    max.overlaps = Inf
    
  ) +
  
  scale_color_manual(
    
    values = c(
      
      "Upregulated" = "#D55E00",
      "Downregulated" = "#0072B2",
      "Not Significant" = "grey70"
      
    )
    
  ) +
  
  labs(
    
    title = "Differentially Expressed Genes",
    
    x = "log2 Fold Change",
    
    y = "-log10 Adjusted P-value",
    
    color = "Expression"
    
  ) +
  
  theme_classic() +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    legend.title = element_text(
      face = "bold"
    )
    
  )


# ----------------------------------------------------------
# 03.6 Save volcano plot
# ----------------------------------------------------------

ggsave(
  
  filename = file.path(
    project_dir,
    "07_figures",
    "Figure1_Volcano.tiff"
  ),
  
  plot = volcano_plot,
  
  width = 6,
  
  height = 6,
  
  dpi = 600,
  
  compression = "lzw"
  
)


# ----------------------------------------------------------
# 03.7 Display volcano plot
# ----------------------------------------------------------

print(volcano_plot)


# ----------------------------------------------------------
# 03.8 Prepare variance-stabilized expression matrix
# ----------------------------------------------------------

expr <- assay(vsd)


# ----------------------------------------------------------
# 03.9 Create display names
# ----------------------------------------------------------

resultsfilter <- resultsfilter %>%
  
  mutate(
    
    Display_Name = ifelse(
      
      is.na(Gene_Symbol) | Gene_Symbol == "",
      
      Ensembl_ID,
      
      Gene_Symbol
      
    )
    
  )


# ----------------------------------------------------------
# 03.10 Calculate composite ranking score
# ----------------------------------------------------------

resultsfilter <- resultsfilter %>%
  
  mutate(
    
    Score =
      abs(log2FoldChange) *
      -log10(padj)
    
  ) %>%
  
  arrange(
    desc(Score)
  )


# ----------------------------------------------------------
# 03.11 Select top 50 DEGs
# ----------------------------------------------------------

Top50_DEGs <- resultsfilter %>%
  
  slice_head(
    n = 50
  )


# ----------------------------------------------------------
# 03.12 Save top 50 DEG table
# ----------------------------------------------------------

write_csv(
  
  Top50_DEGs,
  
  file.path(
    project_dir,
    "10_tables",
    "Top50_DEGs_for_Heatmap.csv"
  )
  
)


# ----------------------------------------------------------
# 03.13 Build top 50 DEG expression matrix
# ----------------------------------------------------------

heatmap_matrix <- expr[
  
  Top50_DEGs$Ensembl_ID,
  
]


# ----------------------------------------------------------
# 03.14 Replace row names with display names
# ----------------------------------------------------------

rownames(heatmap_matrix) <-
  Top50_DEGs$Display_Name


# ----------------------------------------------------------
# 03.15 Rename samples for public repository
# ----------------------------------------------------------

# Original sample accession identifiers are omitted from
# this public repository.

colnames(heatmap_matrix) <- c(
  
  "Case_1",
  "Case_2",
  "Case_3",
  "Control_1",
  "Control_2",
  "Control_3"
  
)


# ----------------------------------------------------------
# 03.16 Z-score normalization
# ----------------------------------------------------------

heatmap_matrix <-
  
  t(
    
    scale(
      t(heatmap_matrix)
    )
    
  )

heatmap_matrix <-
  
  heatmap_matrix[
    
    complete.cases(
      heatmap_matrix
    ),
    
  ]


# ----------------------------------------------------------
# 03.17 Define sample annotation
# ----------------------------------------------------------

ha <- HeatmapAnnotation(
  
  Group = c(
    
    "Case",
    "Case",
    "Case",
    "Control",
    "Control",
    "Control"
    
  ),
  
  col = list(
    
    Group = c(
      
      Case = "#C62828",
      Control = "#1565C0"
      
    )
    
  )
  
)


# ----------------------------------------------------------
# 03.18 Define heatmap color scale
# ----------------------------------------------------------

col_fun <- colorRamp2(
  
  c(-2, 0, 2),
  
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
  
)


# ----------------------------------------------------------
# 03.19 Generate top 50 DEG heatmap
# ----------------------------------------------------------

ht <- Heatmap(
  
  heatmap_matrix,
  
  name = "Z-score",
  
  col = col_fun,
  
  top_annotation = ha,
  
  cluster_rows = TRUE,
  
  cluster_columns = TRUE,
  
  show_row_names = TRUE,
  
  show_column_names = TRUE,
  
  row_names_gp = gpar(
    fontsize = 7
  ),
  
  column_names_gp = gpar(
    fontsize = 10,
    fontface = "bold"
  ),
  
  column_title = "PCOS vs Control",
  
  row_title =
    "Top 50 Differentially Expressed Genes",
  
  heatmap_legend_param = list(
    title = "Expression"
  )
  
)


# ----------------------------------------------------------
# 03.20 Save top 50 DEG heatmap as TIFF
# ----------------------------------------------------------

tiff(
  
  file.path(
    project_dir,
    "07_figures",
    "Figure2_Top50_Heatmap.tiff"
  ),
  
  width = 7,
  
  height = 8,
  
  units = "in",
  
  res = 600,
  
  compression = "lzw"
  
)

draw(ht)

dev.off()


# ----------------------------------------------------------
# 03.21 Save top 50 DEG heatmap as PDF
# ----------------------------------------------------------

pdf(
  
  file.path(
    project_dir,
    "07_figures",
    "Figure2_Top50_Heatmap.pdf"
  ),
  
  width = 7,
  
  height = 8
  
)

draw(ht)

dev.off()


# ----------------------------------------------------------
# 03.22 Prepare lncRNA expression matrix
# ----------------------------------------------------------

lnc_matrix <- expr[
  
  lncRNAs$Ensembl_ID,
  
]


# ----------------------------------------------------------
# 03.23 Rename lncRNA rows
# ----------------------------------------------------------

rownames(lnc_matrix) <-
  lncRNAs$Display_Name


# ----------------------------------------------------------
# 03.24 Rename samples
# ----------------------------------------------------------

# Original sample accession identifiers are omitted from
# this public repository.

colnames(lnc_matrix) <- c(
  
  "Case_1",
  "Case_2",
  "Case_3",
  "Control_1",
  "Control_2",
  "Control_3"
  
)


# ----------------------------------------------------------
# 03.25 Z-score normalization of lncRNA expression
# ----------------------------------------------------------

lnc_matrix <- t(
  
  scale(
    t(lnc_matrix)
  )
  
)

lnc_matrix <-
  
  lnc_matrix[
    
    complete.cases(
      lnc_matrix
    ),
    
  ]


# ----------------------------------------------------------
# 03.26 Generate lncRNA heatmap
# ----------------------------------------------------------

ht_lnc <- Heatmap(
  
  lnc_matrix,
  
  name = "Z-score",
  
  col = col_fun,
  
  top_annotation = ha,
  
  cluster_rows = TRUE,
  
  cluster_columns = TRUE,
  
  show_row_names = TRUE,
  
  show_column_names = TRUE,
  
  row_names_gp = gpar(
    fontsize = 8,
    fontface = "italic"
  ),
  
  column_names_gp = gpar(
    fontsize = 10,
    fontface = "bold"
  ),
  
  row_title = "Differentially Expressed lncRNAs",
  
  column_title = "PCOS vs Control"
  
)


# ----------------------------------------------------------
# 03.27 Save lncRNA heatmap as TIFF
# ----------------------------------------------------------

tiff(
  
  file.path(
    project_dir,
    "07_figures",
    "Figure2B_lncRNA_Heatmap.tiff"
  ),
  
  width = 6,
  
  height = 7,
  
  units = "in",
  
  res = 600,
  
  compression = "lzw"
  
)

draw(ht_lnc)

dev.off()


# ----------------------------------------------------------
# 03.28 Save lncRNA heatmap as PDF
# ----------------------------------------------------------

pdf(
  
  file.path(
    project_dir,
    "07_figures",
    "Figure2B_lncRNA_Heatmap.pdf"
  ),
  
  width = 6,
  
  height = 7
  
)

draw(ht_lnc)

dev.off()
