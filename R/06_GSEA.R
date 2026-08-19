############################################################
# 06. Gene Set Enrichment Analysis (GSEA)
############################################################


# ----------------------------------------------------------
# 06.1 Load required packages
# ----------------------------------------------------------

library(clusterProfiler)
library(ReactomePA)
library(msigdbr)
library(enrichplot)
library(org.Hs.eg.db)
library(readr)
library(dplyr)
library(writexl)


# ----------------------------------------------------------
# 06.2 Create GSEA output directories
# ----------------------------------------------------------

gsea_dir <- file.path(
  project_dir,
  "11_GSEA"
)

dir.create(
  gsea_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(
    gsea_dir,
    "Figures"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(
    gsea_dir,
    "Supplementary"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)


# ----------------------------------------------------------
# 06.3 Prepare ranked gene list
# ----------------------------------------------------------

gsea_df <- results_annotation %>%
  
  filter(
    !is.na(stat),
    !is.na(Entrez_ID)
  ) %>%
  
  distinct(
    Entrez_ID,
    .keep_all = TRUE
  ) %>%
  
  arrange(
    desc(stat)
  )


geneList <- gsea_df$stat

names(geneList) <- gsea_df$Entrez_ID

geneList <- sort(
  geneList,
  decreasing = TRUE
)


# ----------------------------------------------------------
# 06.4 Export ranked gene list
# ----------------------------------------------------------

write_csv(
  
  gsea_df,
  
  file.path(
    
    gsea_dir,
    
    "Ranked_Genes.csv"
    
  )
)


# ----------------------------------------------------------
# 06.5 KEGG GSEA
# ----------------------------------------------------------

gsea_kegg <- gseKEGG(
  
  geneList = geneList,
  
  organism = "hsa",
  
  minGSSize = 10,
  
  maxGSSize = 500,
  
  pvalueCutoff = 0.05,
  
  verbose = FALSE
  
)


# ----------------------------------------------------------
# 06.6 Reactome GSEA
# ----------------------------------------------------------

gsea_reactome <- gsePathway(
  
  geneList = geneList,
  
  organism = "human",
  
  minGSSize = 10,
  
  maxGSSize = 500,
  
  pvalueCutoff = 0.05,
  
  verbose = FALSE
  
)


# ----------------------------------------------------------
# 06.7 Hallmark GSEA
# ----------------------------------------------------------

hallmark_t2g <-
  
  msigdbr(
    
    species = "Homo sapiens",
    
    collection = "H"
    
  ) %>%
  
  dplyr::select(
    
    gs_name,
    
    ncbi_gene
    
  ) %>%
  
  dplyr::filter(
    
    !is.na(ncbi_gene)
    
  )


gsea_hallmark <- GSEA(
  
  geneList,
  
  TERM2GENE = hallmark_t2g,
  
  minGSSize = 10,
  
  maxGSSize = 500,
  
  pvalueCutoff = 0.05,
  
  verbose = FALSE
  
)


# ----------------------------------------------------------
# 06.8 Prepare Hallmark GSEA result table
# ----------------------------------------------------------

Hallmark_GSEA_Table <-
  
  as.data.frame(
    gsea_hallmark
  ) %>%
  
  arrange(
    
    p.adjust,
    
    desc(abs(NES))
    
  )


# ----------------------------------------------------------
# 06.9 Inspect Hallmark GSEA results
# ----------------------------------------------------------

head(
  hallmark_t2g
)

Hallmark_GSEA_Table %>%
  
  dplyr::select(
    
    Description,
    NES,
    p.adjust,
    setSize
    
  ) %>%
  
  head(20)


hallmark_t2g$ncbi_gene <-
  as.character(
    hallmark_t2g$ncbi_gene
  )


head(
  Hallmark_GSEA_Table[
    
    ,
    
    c(
      "Description",
      "NES",
      "p.adjust",
      "setSize"
    )
    
  ],
  
  20
  
)


# ----------------------------------------------------------
# 06.10 Save GSEA objects
# ----------------------------------------------------------

saveRDS(
  
  gsea_kegg,
  
  file.path(
    
    gsea_dir,
    
    "GSEA_KEGG.rds"
    
  )
)


saveRDS(
  
  gsea_reactome,
  
  file.path(
    
    gsea_dir,
    
    "GSEA_Reactome.rds"
    
  )
)


saveRDS(
  
  gsea_hallmark,
  
  file.path(
    
    gsea_dir,
    
    "GSEA_Hallmark.rds"
    
  )
)


# ----------------------------------------------------------
# 06.11 Prepare GSEA result tables
# ----------------------------------------------------------

KEGG_GSEA_Table <-
  
  as.data.frame(
    gsea_kegg
  ) %>%
  
  arrange(
    
    p.adjust,
    
    desc(abs(NES))
    
  )


Reactome_GSEA_Table <-
  
  as.data.frame(
    gsea_reactome
  ) %>%
  
  arrange(
    
    p.adjust,
    
    desc(abs(NES))
    
  )


Hallmark_GSEA_Table <-
  
  as.data.frame(
    gsea_hallmark
  ) %>%
  
  arrange(
    
    p.adjust,
    
    desc(abs(NES))
    
  )


# ----------------------------------------------------------
# 06.12 Export supplementary GSEA tables
# ----------------------------------------------------------

write_xlsx(
  
  KEGG_GSEA_Table,
  
  file.path(
    
    gsea_dir,
    
    "Supplementary",
    
    "Supplementary_Table_S11_GSEA_KEGG.xlsx"
    
  )
)


write_xlsx(
  
  Reactome_GSEA_Table,
  
  file.path(
    
    gsea_dir,
    
    "Supplementary",
    
    "Supplementary_Table_S12_GSEA_Reactome.xlsx"
    
  )
)


write_xlsx(
  
  Hallmark_GSEA_Table,
  
  file.path(
    
    gsea_dir,
    
    "Supplementary",
    
    "Supplementary_Table_S13_GSEA_Hallmark.xlsx"
    
  )
)


# ----------------------------------------------------------
# 06.13 Prepare Hallmark figure data
# ----------------------------------------------------------

keep_hallmark <- c(
  
  "HALLMARK_INFLAMMATORY_RESPONSE",
  
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  
  "HALLMARK_COMPLEMENT",
  
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  
  "HALLMARK_MYOGENESIS",
  
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  
  "HALLMARK_FATTY_ACID_METABOLISM",
  
  "HALLMARK_MTORC1_SIGNALING"
  
)


hallmark_plot <-
  
  Hallmark_GSEA_Table %>%
  
  filter(
    
    Description %in% keep_hallmark
    
  ) %>%
  
  mutate(
    
    Direction = ifelse(
      
      NES > 0,
      
      "Activated",
      
      "Suppressed"
      
    ),
    
    FDR = -log10(
      p.adjust
    ),
    
    Description = gsub(
      
      "HALLMARK_",
      
      "",
      
      Description
      
    ),
    
    Description = gsub(
      
      "_",
      
      " ",
      
      Description
      
    ),
    
    Description = stringr::str_to_title(
      Description
    ),
    
    Description = stringr::str_wrap(
      
      Description,
      
      35
      
    )
    
  )


# ----------------------------------------------------------
# 06.14 Generate Hallmark NES plot
# ----------------------------------------------------------

Fig4A <-
  
  ggplot(
    
    hallmark_plot,
    
    aes(
      
      x = NES,
      
      y = reorder(
        Description,
        NES
      )
      
    )
    
  ) +
  
  geom_vline(
    
    xintercept = 0,
    
    linetype = 2,
    
    colour = "grey60"
    
  ) +
  
  geom_segment(
    
    aes(
      
      x = 0,
      
      xend = NES,
      
      y = Description,
      
      yend = Description
      
    ),
    
    linewidth = 0.8,
    
    colour = "grey75"
    
  ) +
  
  geom_point(
    
    aes(
      
      size = FDR,
      
      fill = Direction
      
    ),
    
    shape = 21,
    
    colour = "black",
    
    stroke = 0.3
    
  ) +
  
  scale_fill_manual(
    
    values = c(
      
      "Activated" = "#B40426",
      
      "Suppressed" = "#3B4CC0"
      
    ),
    
    name = NULL
    
  ) +
  
  scale_size(
    
    range = c(
      4,
      10
    ),
    
    name = expression(
      -log[10](FDR)
    )
    
  ) +
  
  labs(
    
    x = "Normalized Enrichment Score (NES)",
    
    y = NULL
    
  ) +
  
  theme_classic(
    
    base_family = "Arial",
    
    base_size = 11
    
  ) +
  
  theme(
    
    axis.text = element_text(
      colour = "black"
    ),
    
    legend.title = element_text(
      face = "bold"
    )
    
  )


# ----------------------------------------------------------
# 06.15 Generate Hallmark enrichment plots
# ----------------------------------------------------------

Fig4B <-
  
  gseaplot2(
    
    gsea_hallmark,
    
    geneSetID =
      "HALLMARK_INFLAMMATORY_RESPONSE",
    
    title =
      "Inflammatory Response",
    
    pvalue_table = TRUE,
    
    base_size = 11
    
  )


Fig4C <-
  
  gseaplot2(
    
    gsea_hallmark,
    
    geneSetID =
      "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
    
    title =
      "Cholesterol Homeostasis",
    
    pvalue_table = TRUE,
    
    base_size = 11
    
  )


# ----------------------------------------------------------
# 06.16 Define GSEA figure directory
# ----------------------------------------------------------

figure_dir <- file.path(
  
  gsea_dir,
  
  "Figures"
  
)


# ----------------------------------------------------------
# 06.17 Save Hallmark NES figure
# ----------------------------------------------------------

ggsave(
  
  file.path(
    
    figure_dir,
    
    "Fig4A_Hallmark_NES.tiff"
    
  ),
  
  Fig4A,
  
  dpi = 600,
  
  width = 7,
  
  height = 5.5,
  
  compression = "lzw"
  
)


# ----------------------------------------------------------
# 06.18 Save inflammatory response GSEA figure
# ----------------------------------------------------------

ggsave(
  
  file.path(
    
    figure_dir,
    
    "Fig4B_Inflammatory_Response.tiff"
    
  ),
  
  Fig4B,
  
  dpi = 600,
  
  width = 15,
  
  height = 5,
  
  compression = "lzw"
  
)


# ----------------------------------------------------------
# 06.19 Save cholesterol homeostasis GSEA figure
# ----------------------------------------------------------

ggsave(
  
  file.path(
    
    figure_dir,
    
    "Fig4C_Cholesterol_Homeostasis.tiff"
    
  ),
  
  Fig4C,
  
  dpi = 600,
  
  width = 15,
  
  height = 5,
  
  compression = "lzw"
  
)


# ----------------------------------------------------------
# 06.20 Inspect gene biotypes
# ----------------------------------------------------------

table(
  resultsfilter$Gene_Biotype
)


resultsfilter %>%
  
  dplyr::count(
    
    Gene_Biotype,
    
    sort = TRUE
    
  )
