############################################################
# 09. lncRNA-Hub Gene Correlation Analysis
############################################################


# ----------------------------------------------------------
# 09.1 Load required packages
# ----------------------------------------------------------

library(dplyr)
library(Hmisc)
library(readr)
library(tibble)
library(tidyr)
library(ComplexHeatmap)
library(circlize)


# ----------------------------------------------------------
# 09.2 Create output directories
# ----------------------------------------------------------

correlation_dir <- file.path(
  project_dir,
  "13_correlation network_analysis"
)

dir.create(
  correlation_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  file.path(
    correlation_dir,
    "Correlation"
  ),
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  file.path(
    correlation_dir,
    "Figures"
  ),
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  file.path(
    correlation_dir,
    "Tables"
  ),
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  file.path(
    correlation_dir,
    "Network_for_Cytoscape"
  ),
  showWarnings = FALSE,
  recursive = TRUE
)


# ----------------------------------------------------------
# 09.3 Load normalized expression matrix
# ----------------------------------------------------------

# The variance stabilizing transformation (VST) was
# already performed in the previous section.
# The existing normalized expression matrix is therefore
# reused to avoid redundant computation.

expr_matrix <- assay(vsd)


# ----------------------------------------------------------
# 09.4 Prepare consensus hub gene list
# ----------------------------------------------------------

hub_genes <- Consensus_Hubs$Gene_Symbol


# ----------------------------------------------------------
# 09.5 Prepare DE lncRNA list
# ----------------------------------------------------------

lncRNA_ids <- lncRNAs$Ensembl_ID


# ----------------------------------------------------------
# 09.6 Summarize input features
# ----------------------------------------------------------

cat(
  
  "Consensus hub genes:",
  
  length(hub_genes),
  
  "\n"
  
)

cat(
  
  "Differentially expressed lncRNAs:",
  
  length(lncRNA_ids),
  
  "\n"
  
)


# ----------------------------------------------------------
# 09.7 Convert expression matrix to data frame
# ----------------------------------------------------------

expr_df <- as.data.frame(
  expr_matrix
)

expr_df$Ensembl_ID <- rownames(
  expr_df
)


# ----------------------------------------------------------
# 09.8 Add DEG annotation
# ----------------------------------------------------------

# resultsfilter contains the significant DEG set,
# including lncRNAs without official gene symbols.

annotation_df <- resultsfilter %>%
  
  dplyr::select(
    
    Ensembl_ID,
    
    Display_Name,
    
    Gene_Symbol,
    
    Gene_Biotype
    
  )


expr_df <- expr_df %>%
  
  left_join(
    
    annotation_df,
    
    by = "Ensembl_ID"
    
  )


# ----------------------------------------------------------
# 09.9 Keep annotated genes
# ----------------------------------------------------------

expr_df <- expr_df %>%
  
  filter(
    
    !is.na(Display_Name)
    
  )


# ----------------------------------------------------------
# 09.10 Extract consensus hub-gene expression
# ----------------------------------------------------------

hub_expression <- expr_df %>%
  
  filter(
    
    Gene_Symbol %in% hub_genes
    
  )


# ----------------------------------------------------------
# 09.11 Extract DE lncRNA expression
# ----------------------------------------------------------

lnc_expression <- expr_df %>%
  
  filter(
    
    Ensembl_ID %in% lncRNAs$Ensembl_ID
    
  )


# ----------------------------------------------------------
# 09.12 Check extracted features
# ----------------------------------------------------------

cat(
  
  "\nConsensus hub genes found:",
  
  nrow(hub_expression),
  
  "\n"
  
)

cat(
  
  "DE lncRNAs found:",
  
  nrow(lnc_expression),
  
  "\n"
  
)


# ----------------------------------------------------------
# 09.13 Export expression matrices
# ----------------------------------------------------------

write_csv(
  
  hub_expression,
  
  file.path(
    
    correlation_dir,
    
    "Correlation",
    
    "Hub_Gene_Expression.csv"
    
  )
)


write_csv(
  
  lnc_expression,
  
  file.path(
    
    correlation_dir,
    
    "Correlation",
    
    "lncRNA_Expression.csv"
    
  )
)


# ----------------------------------------------------------
# 09.14 Prepare matrices for correlation analysis
# ----------------------------------------------------------

hub_matrix <- hub_expression %>%
  
  column_to_rownames(
    "Gene_Symbol"
  ) %>%
  
  dplyr::select(
    
    -Ensembl_ID,
    
    -Display_Name,
    
    -Gene_Biotype
    
  )


lnc_matrix <- lnc_expression %>%
  
  column_to_rownames(
    "Display_Name"
  ) %>%
  
  dplyr::select(
    
    -Ensembl_ID,
    
    -Gene_Symbol,
    
    -Gene_Biotype
    
  )


# ----------------------------------------------------------
# 09.15 Check matrix dimensions
# ----------------------------------------------------------

cat(
  
  "\nHub matrix:",
  
  dim(hub_matrix),
  
  "\n"
  
)

cat(
  
  "lncRNA matrix:",
  
  dim(lnc_matrix),
  
  "\n"
  
)


# Optional consistency check

sum(
  rownames(expr_matrix) %in%
    resultsfilter$Ensembl_ID
)

head(
  rownames(expr_matrix)
)


# ----------------------------------------------------------
# 09.16 Initialize correlation result table
# ----------------------------------------------------------

Correlation_Results <- data.frame()


# ----------------------------------------------------------
# 09.17 Calculate Spearman correlations
# ----------------------------------------------------------

for (i in rownames(lnc_matrix)) {
  
  for (j in rownames(hub_matrix)) {
    
    cor_result <- rcorr(
      
      as.numeric(
        lnc_matrix[i, ]
      ),
      
      as.numeric(
        hub_matrix[j, ]
      ),
      
      type = "spearman"
      
    )
    
    
    Correlation_Results <- rbind(
      
      Correlation_Results,
      
      data.frame(
        
        lncRNA = i,
        
        Hub_Gene = j,
        
        Spearman_Rho =
          cor_result$r[1, 2],
        
        P_value =
          cor_result$P[1, 2]
        
      )
      
    )
    
  }
  
}


# ----------------------------------------------------------
# 09.18 Adjust P-values using Benjamini-Hochberg
# ----------------------------------------------------------

Correlation_Results$FDR <-
  
  p.adjust(
    
    Correlation_Results$P_value,
    
    method = "BH"
    
  )


# ----------------------------------------------------------
# 09.19 Determine correlation direction
# ----------------------------------------------------------

Correlation_Results <- Correlation_Results %>%
  
  mutate(
    
    Direction = case_when(
      
      Spearman_Rho > 0 ~ "Positive",
      
      Spearman_Rho < 0 ~ "Negative",
      
      TRUE ~ "None"
      
    )
    
  )


# ----------------------------------------------------------
# 09.20 Rank correlations
# ----------------------------------------------------------

Correlation_Results <- Correlation_Results %>%
  
  arrange(
    
    FDR,
    
    desc(
      abs(Spearman_Rho)
    )
    
  )


# ----------------------------------------------------------
# 09.21 Save complete correlation table
# ----------------------------------------------------------

write_csv(
  
  Correlation_Results,
  
  file.path(
    
    correlation_dir,
    
    "Correlation",
    
    "All_lncRNA_HubGene_Correlations.csv"
    
  )
)


# ----------------------------------------------------------
# 09.22 Summarize correlation results
# ----------------------------------------------------------

cat(
  
  "\nTotal correlations:",
  
  nrow(Correlation_Results),
  
  "\n"
  
)

cat(
  
  "Unique lncRNAs:",
  
  length(
    unique(
      Correlation_Results$lncRNA
    )
  ),
  
  "\n"
  
)

cat(
  
  "Unique hub genes:",
  
  length(
    unique(
      Correlation_Results$Hub_Gene
    )
  ),
  
  "\n"
  
)

summary(
  Correlation_Results$Spearman_Rho
)

summary(
  Correlation_Results$P_value
)

summary(
  Correlation_Results$FDR
)


# ----------------------------------------------------------
# 09.23 Inspect strong correlations
# ----------------------------------------------------------

Correlation_Results %>%
  
  filter(
    
    abs(Spearman_Rho) >= 0.8,
    
    P_value < 0.05
    
  ) %>%
  
  nrow()


Correlation_Results %>%
  
  filter(
    
    abs(Spearman_Rho) >= 0.8,
    
    P_value < 0.05
    
  )


# ----------------------------------------------------------
# 09.24 Filter significant correlations
# ----------------------------------------------------------

Correlation_Significant <-
  
  Correlation_Results %>%
  
  filter(
    
    abs(Spearman_Rho) >= 0.8,
    
    P_value < 0.05
    
  ) %>%
  
  arrange(
    
    desc(
      abs(Spearman_Rho)
    )
    
  )


cat(
  
  "\nSignificant correlations:",
  
  nrow(
    Correlation_Significant
  ),
  
  "\n"
  
)


# ----------------------------------------------------------
# 09.25 Save significant correlations
# ----------------------------------------------------------

write_csv(
  
  Correlation_Significant,
  
  file.path(
    
    correlation_dir,
    
    "Correlation",
    
    "Significant_Hub_lncRNA_Correlations.csv"
    
  )
)


# ----------------------------------------------------------
# 09.26 Save all correlations
# ----------------------------------------------------------

write_csv(
  
  Correlation_Results,
  
  file.path(
    
    correlation_dir,
    
    "Correlation",
    
    "All_Hub_lncRNA_Correlations.csv"
    
  )
)


# ----------------------------------------------------------
# 09.27 Prepare Cytoscape edge table
# ----------------------------------------------------------

Edge_Table <- Correlation_Significant %>%
  
  transmute(
    
    Source = lncRNA,
    
    Target = Hub_Gene,
    
    Weight = abs(
      Spearman_Rho
    ),
    
    Correlation = Spearman_Rho,
    
    Direction = Direction,
    
    P_value = P_value
    
  )


# ----------------------------------------------------------
# 09.28 Export Cytoscape edge table
# ----------------------------------------------------------

write_csv(
  
  Edge_Table,
  
  file.path(
    
    correlation_dir,
    
    "Network_for_Cytoscape",
    
    "Cytoscape_Edges.csv"
    
  )
)


cat(
  
  "\nEdges exported:",
  
  nrow(
    Edge_Table
  ),
  
  "\n"
  
)


# ----------------------------------------------------------
# 09.29 Prepare Cytoscape node table
# ----------------------------------------------------------

Hub_nodes <-
  
  tibble(
    
    Node = unique(
      Edge_Table$Target
    ),
    
    Type = "Hub_gene"
    
  )


lnc_nodes <-
  
  tibble(
    
    Node = unique(
      Edge_Table$Source
    ),
    
    Type = "lncRNA"
    
  )


Node_Table <-
  
  bind_rows(
    
    Hub_nodes,
    
    lnc_nodes
    
  )


# ----------------------------------------------------------
# 09.30 Export Cytoscape node table
# ----------------------------------------------------------

write_csv(
  
  Node_Table,
  
  file.path(
    
    correlation_dir,
    
    "Network_for_Cytoscape",
    
    "Cytoscape_Nodes.csv"
    
  )
)


cat(
  
  "\nNodes exported:",
  
  nrow(
    Node_Table
  ),
  
  "\n"
  
)


# ----------------------------------------------------------
# 09.31 Generate correlation summary
# ----------------------------------------------------------

Correlation_Summary <-
  
  Correlation_Significant %>%
  
  group_by(
    Direction
  ) %>%
  
  summarise(
    
    Number_of_Correlations = n(),
    
    Mean_Rho = mean(
      Spearman_Rho
    ),
    
    Median_Rho = median(
      Spearman_Rho
    ),
    
    .groups = "drop"
    
  )


print(
  Correlation_Summary
)


# ----------------------------------------------------------
# 09.32 Export correlation summary
# ----------------------------------------------------------

write_csv(
  
  Correlation_Summary,
  
  file.path(
    
    correlation_dir,
    
    "Correlation",
    
    "Correlation_Summary.csv"
    
  )
)


# ----------------------------------------------------------
# 09.33 Prepare correlation heatmap data
# ----------------------------------------------------------

Correlation_Filtered <-
  
  Correlation_Results %>%
  
  filter(
    
    abs(Spearman_Rho) >= 0.80,
    
    P_value < 0.05
    
  )


# ----------------------------------------------------------
# 09.34 Create correlation matrix
# ----------------------------------------------------------

Correlation_Matrix <-
  
  Correlation_Filtered %>%
  
  dplyr::select(
    
    lncRNA,
    
    Hub_Gene,
    
    Spearman_Rho
    
  ) %>%
  
  pivot_wider(
    
    names_from = Hub_Gene,
    
    values_from = Spearman_Rho
    
  )


# ----------------------------------------------------------
# 09.35 Convert correlation matrix to matrix object
# ----------------------------------------------------------

rho_matrix <-
  
  as.data.frame(
    Correlation_Matrix
  )


rownames(rho_matrix) <-
  rho_matrix$lncRNA

rho_matrix$lncRNA <- NULL

rho_matrix <-
  
  as.matrix(
    rho_matrix
  )


# ----------------------------------------------------------
# 09.36 Replace lncRNA Ensembl IDs with display names
# ----------------------------------------------------------

lnc_annotation <-
  
  lncRNAs %>%
  
  dplyr::select(
    
    Ensembl_ID,
    
    Display_Name
    
  )


rownames(rho_matrix) <-
  
  sapply(
    
    rownames(rho_matrix),
    
    function(x) {
      
      if (
        x %in%
        lnc_annotation$Ensembl_ID
      ) {
        
        lnc_annotation$Display_Name[
          
          match(
            
            x,
            
            lnc_annotation$Ensembl_ID
            
          )
          
        ]
        
      } else {
        
        x
        
      }
      
    }
    
  )


# ----------------------------------------------------------
# 09.37 Remove empty rows and columns
# ----------------------------------------------------------

rho_matrix <-
  
  rho_matrix[
    
    rowSums(
      !is.na(rho_matrix)
    ) > 0,
    
  ]


rho_matrix <-
  
  rho_matrix[
    
    ,
    
    colSums(
      !is.na(rho_matrix)
    ) > 0
    
  ]


# ----------------------------------------------------------
# 09.38 Calculate node degree
# ----------------------------------------------------------

row_degree <-
  
  rowSums(
    !is.na(rho_matrix)
  )


col_degree <-
  
  colSums(
    !is.na(rho_matrix)
  )


# ----------------------------------------------------------
# 09.39 Order heatmap rows and columns
# ----------------------------------------------------------

rho_matrix <-
  
  rho_matrix[
    
    order(
      
      row_degree,
      
      decreasing = TRUE
      
    ),
    
  ]


rho_matrix <-
  
  rho_matrix[
    
    ,
    
    order(
      
      col_degree,
      
      decreasing = TRUE
      
    )
    
  ]


# ----------------------------------------------------------
# 09.40 Define correlation color scale
# ----------------------------------------------------------

col_fun <-
  
  colorRamp2(
    
    c(
      -1,
      0,
      1
    ),
    
    c(
      
      "#2166AC",
      
      "white",
      
      "#B2182B"
      
    )
    
  )


# ----------------------------------------------------------
# 09.41 Generate correlation heatmap
# ----------------------------------------------------------

ht_corr <-
  
  Heatmap(
    
    rho_matrix,
    
    name = "Spearman rho",
    
    col = col_fun,
    
    na_col = "white",
    
    border = TRUE,
    
    rect_gp = gpar(
      
      col = "grey85",
      
      lwd = 0.3
      
    ),
    
    cluster_rows = FALSE,
    
    cluster_columns = FALSE,
    
    show_row_names = TRUE,
    
    show_column_names = TRUE,
    
    row_names_gp = gpar(
      
      fontsize = 10,
      
      fontface = "italic"
      
    ),
    
    column_names_gp = gpar(
      
      fontsize = 10,
      
      fontface = "bold"
      
    ),
    
    column_names_rot = 45,
    
    row_title =
      "Differentially Expressed lncRNAs",
    
    column_title =
      "Hub Genes",
    
    heatmap_legend_param = list(
      
      title = expression(rho),
      
      at = c(
        -1,
        -0.5,
        0,
        0.5,
        1
      )
      
    )
    
  )


# ----------------------------------------------------------
# 09.42 Save correlation heatmap as TIFF
# ----------------------------------------------------------

tiff(
  
  file.path(
    
    correlation_dir,
    
    "Figures",
    
    "Figure6A_Correlation_Heatmap.tiff"
    
  ),
  
  width = 8,
  
  height = 7,
  
  units = "in",
  
  res = 600,
  
  compression = "lzw"
  
)

draw(
  ht_corr
)

dev.off()


# ----------------------------------------------------------
# 09.43 Save correlation heatmap as PDF
# ----------------------------------------------------------

pdf(
  
  file.path(
    
    correlation_dir,
    
    "Figures",
    
    "Figure6A_Correlation_Heatmap.pdf"
    
  ),
  
  width = 8,
  
  height = 7
  
)

draw(
  ht_corr
)

dev.off()
