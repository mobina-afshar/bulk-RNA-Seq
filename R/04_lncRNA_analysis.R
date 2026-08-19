############################################################
# 04. Differentially Expressed lncRNA Analysis
############################################################


# ----------------------------------------------------------
# 04.1 Define lncRNA biotypes
# ----------------------------------------------------------

lncRNA_biotypes <- c(
  
  "lncRNA",
  "antisense",
  "sense_intronic",
  "sense_overlapping",
  "processed_transcript",
  "3prime_overlapping_ncRNA",
  "bidirectional_promoter_lncRNA",
  "macro_lncRNA"
  
)


# ----------------------------------------------------------
# 04.2 Extract differentially expressed lncRNAs
# ----------------------------------------------------------

lncRNAs <- resultsfilter %>%
  
  dplyr::filter(
    
    Gene_Biotype %in% lncRNA_biotypes
    
  )


# ----------------------------------------------------------
# 04.3 Create display names
# ----------------------------------------------------------

lncRNAs <- lncRNAs %>%
  
  dplyr::mutate(
    
    Display_Name = dplyr::if_else(
      
      is.na(Gene_Symbol) | Gene_Symbol == "",
      
      Ensembl_ID,
      
      Gene_Symbol
      
    )
    
  )


# ----------------------------------------------------------
# 04.4 Add expression direction
# ----------------------------------------------------------

lncRNAs <- lncRNAs %>%
  
  dplyr::mutate(
    
    Direction = dplyr::if_else(
      
      log2FoldChange > 0,
      
      "Upregulated",
      
      "Downregulated"
      
    )
    
  )


# ----------------------------------------------------------
# 04.5 Sort lncRNAs by statistical significance
# ----------------------------------------------------------

lncRNAs <- lncRNAs %>%
  
  dplyr::arrange(
    padj
  )


# ----------------------------------------------------------
# 04.6 Separate upregulated and downregulated lncRNAs
# ----------------------------------------------------------

up_lncRNAs <- lncRNAs %>%
  
  dplyr::filter(
    Direction == "Upregulated"
  )


down_lncRNAs <- lncRNAs %>%
  
  dplyr::filter(
    Direction == "Downregulated"
  )


# ----------------------------------------------------------
# 04.7 Generate summary statistics
# ----------------------------------------------------------

cat(
  
  "\nTotal DE lncRNAs:",
  
  nrow(lncRNAs),
  
  "\n"
  
)

cat(
  
  "Upregulated lncRNAs:",
  
  nrow(up_lncRNAs),
  
  "\n"
  
)

cat(
  
  "Downregulated lncRNAs:",
  
  nrow(down_lncRNAs),
  
  "\n"
  
)


# ----------------------------------------------------------
# 04.8 Summarize lncRNA biotypes
# ----------------------------------------------------------

lncRNA_biotype_summary <- lncRNAs %>%
  
  dplyr::count(
    
    Gene_Biotype,
    
    sort = TRUE
    
  )


# ----------------------------------------------------------
# 04.9 Export lncRNA biotype summary
# ----------------------------------------------------------

write_csv(
  
  lncRNA_biotype_summary,
  
  file.path(
    
    project_dir,
    
    "08_lncRNA",
    
    "lncRNA_Biotype_Summary.csv"
    
  )
  
)


# ----------------------------------------------------------
# 04.10 Export complete lncRNA table
# ----------------------------------------------------------

write_csv(
  
  lncRNAs,
  
  file.path(
    
    project_dir,
    
    "08_lncRNA",
    
    "Differentially_Expressed_lncRNAs.csv"
    
  )
  
)


# ----------------------------------------------------------
# 04.11 Export upregulated lncRNAs
# ----------------------------------------------------------

write_csv(
  
  up_lncRNAs,
  
  file.path(
    
    project_dir,
    
    "08_lncRNA",
    
    "Upregulated_lncRNAs.csv"
    
  )
  
)


# ----------------------------------------------------------
# 04.12 Export downregulated lncRNAs
# ----------------------------------------------------------

write_csv(
  
  down_lncRNAs,
  
  file.path(
    
    project_dir,
    
    "08_lncRNA",
    
    "Downregulated_lncRNAs.csv"
    
  )
  
)


# ----------------------------------------------------------
# 04.13 Save lncRNA R objects
# ----------------------------------------------------------

saveRDS(
  
  lncRNAs,
  
  file.path(
    
    project_dir,
    
    "08_lncRNA",
    
    "Differentially_Expressed_lncRNAs.rds"
    
  )
)


saveRDS(
  
  up_lncRNAs,
  
  file.path(
    
    project_dir,
    
    "08_lncRNA",
    
    "Upregulated_lncRNAs.rds"
    
  )
)


saveRDS(
  
  down_lncRNAs,
  
  file.path(
    
    project_dir,
    
    "08_lncRNA",
    
    "Downregulated_lncRNAs.rds"
    
  )
)


# ----------------------------------------------------------
# 04.14 Prepare supplementary lncRNA table
# ----------------------------------------------------------

Supplementary_lncRNAs <- lncRNAs %>%
  
  dplyr::select(
    
    Ensembl_ID,
    
    Display_Name,
    
    Gene_Symbol,
    
    Gene_Biotype,
    
    Direction,
    
    baseMean,
    
    log2FoldChange,
    
    lfcSE,
    
    stat,
    
    pvalue,
    
    padj
    
  )


# ----------------------------------------------------------
# 04.15 Export supplementary lncRNA table
# ----------------------------------------------------------

write_csv(
  
  Supplementary_lncRNAs,
  
  file.path(
    
    project_dir,
    
    "10_tables",
    
    "Supplementary_Table_DE_lncRNAs.csv"
    
  )
  
)


# ----------------------------------------------------------
# 04.16 Identify top lncRNAs
# ----------------------------------------------------------

Top_lncRNAs <- Supplementary_lncRNAs %>%
  
  dplyr::slice_head(
    n = 20
  )


# ----------------------------------------------------------
# 04.17 Display top lncRNAs
# ----------------------------------------------------------

print(
  Top_lncRNAs
)


# ----------------------------------------------------------
# 04.18 Replace missing gene symbols
# ----------------------------------------------------------

# For lncRNAs without an official gene symbol,
# use the display name for downstream visualization
# and reporting.

lncRNAs$Gene_Symbol[
  
  is.na(lncRNAs$Gene_Symbol) |
    lncRNAs$Gene_Symbol == ""
  
] <-
  
  lncRNAs$Display_Name[
    
    is.na(lncRNAs$Gene_Symbol) |
      lncRNAs$Gene_Symbol == ""
    
  ]
