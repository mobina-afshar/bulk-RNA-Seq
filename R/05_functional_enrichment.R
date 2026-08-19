############################################################
# 05. Functional Enrichment Analysis
############################################################


# ----------------------------------------------------------
# 05.1 Load required packages
# ----------------------------------------------------------

library(enrichR)
library(dplyr)
library(readr)
library(openxlsx)


# ----------------------------------------------------------
# 05.2 Check Enrichr server
# ----------------------------------------------------------

websiteLive <- getOption(
  "enrichR.live"
)

if (!websiteLive) {
  
  stop(
    "Enrichr server is not available."
  )
  
}

listEnrichrSites()

setEnrichrSite(
  "Enrichr"
)


# ----------------------------------------------------------
# 05.3 Create enrichment output directories
# ----------------------------------------------------------

dir.create(
  
  file.path(
    project_dir,
    "09_enrichment",
    "GO"
  ),
  
  recursive = TRUE,
  
  showWarnings = FALSE
  
)

dir.create(
  
  file.path(
    project_dir,
    "09_enrichment",
    "KEGG"
  ),
  
  recursive = TRUE,
  
  showWarnings = FALSE
  
)

dir.create(
  
  file.path(
    project_dir,
    "09_enrichment",
    "Reactome"
  ),
  
  recursive = TRUE,
  
  showWarnings = FALSE
  
)


# ----------------------------------------------------------
# 05.4 Prepare upregulated and downregulated gene lists
# ----------------------------------------------------------

up_genes <- up_regulates %>%
  
  filter(
    
    !is.na(Gene_Symbol),
    
    Gene_Symbol != ""
    
  ) %>%
  
  pull(
    Gene_Symbol
  ) %>%
  
  unique()


down_genes <- down_regulates %>%
  
  filter(
    
    !is.na(Gene_Symbol),
    
    Gene_Symbol != ""
    
  ) %>%
  
  pull(
    Gene_Symbol
  ) %>%
  
  unique()


cat(
  
  "\nUpregulated genes:",
  
  length(up_genes),
  
  "\n"
  
)

cat(
  
  "Downregulated genes:",
  
  length(down_genes),
  
  "\n"
  
)


# ----------------------------------------------------------
# 05.5 Define enrichment databases
# ----------------------------------------------------------

GO_dbs <- c(
  
  "GO_Biological_Process_2026",
  
  "GO_Molecular_Function_2026",
  
  "GO_Cellular_Component_2026"
  
)

KEGG_db <- c(
  
  "KEGG_2026"
  
)

Reactome_db <- c(
  
  "Reactome_Pathways_2024"
  
)


# ----------------------------------------------------------
# 05.6 Run GO enrichment
# ----------------------------------------------------------

GO_up <- enrichr(
  
  up_genes,
  
  GO_dbs
  
)

GO_down <- enrichr(
  
  down_genes,
  
  GO_dbs
  
)


# ----------------------------------------------------------
# 05.7 Run KEGG enrichment
# ----------------------------------------------------------

KEGG_up <- enrichr(
  
  up_genes,
  
  KEGG_db
  
)

KEGG_down <- enrichr(
  
  down_genes,
  
  KEGG_db
  
)


# ----------------------------------------------------------
# 05.8 Run Reactome enrichment
# ----------------------------------------------------------

Reactome_up <- enrichr(
  
  up_genes,
  
  Reactome_db
  
)

Reactome_down <- enrichr(
  
  down_genes,
  
  Reactome_db
  
)


# ----------------------------------------------------------
# 05.9 Save GO tables
# ----------------------------------------------------------

for(i in names(GO_up)){
  
  write_csv(
    
    GO_up[[i]],
    
    file.path(
      
      project_dir,
      
      "09_enrichment",
      
      "GO",
      
      paste0(
        i,
        "_Up.csv"
      )
      
    )
    
  )
  
}


for(i in names(GO_down)){
  
  write_csv(
    
    GO_down[[i]],
    
    file.path(
      
      project_dir,
      
      "09_enrichment",
      
      "GO",
      
      paste0(
        i,
        "_Down.csv"
      )
      
    )
    
  )
  
}


# ----------------------------------------------------------
# 05.10 Save KEGG tables
# ----------------------------------------------------------

write_csv(
  
  KEGG_up[[1]],
  
  file.path(
    
    project_dir,
    
    "09_enrichment",
    
    "KEGG",
    
    "KEGG_Up.csv"
    
  )
  
)

write_csv(
  
  KEGG_down[[1]],
  
  file.path(
    
    project_dir,
    
    "09_enrichment",
    
    "KEGG",
    
    "KEGG_Down.csv"
    
  )
  
)


# ----------------------------------------------------------
# 05.11 Save Reactome tables
# ----------------------------------------------------------

write_csv(
  
  Reactome_up[[1]],
  
  file.path(
    
    project_dir,
    
    "09_enrichment",
    
    "Reactome",
    
    "Reactome_Up.csv"
    
  )
  
)

write_csv(
  
  Reactome_down[[1]],
  
  file.path(
    
    project_dir,
    
    "09_enrichment",
    
    "Reactome_Down.csv"
    
  )
  
)


# ----------------------------------------------------------
# 05.12 Create supplementary Excel workbook
# ----------------------------------------------------------

wb <- createWorkbook()


# ----------------------------------------------------------
# 05.13 Add GO results to supplementary workbook
# ----------------------------------------------------------

for(i in names(GO_up)){
  
  addWorksheet(
    
    wb,
    
    paste0(
      i,
      "_Up"
    )
    
  )
  
  writeData(
    
    wb,
    
    sheet = paste0(
      i,
      "_Up"
    ),
    
    GO_up[[i]]
    
  )
  
}


for(i in names(GO_down)){
  
  addWorksheet(
    
    wb,
    
    paste0(
      i,
      "_Down"
    )
    
  )
  
  writeData(
    
    wb,
    
    sheet = paste0(
      i,
      "_Down"
    ),
    
    GO_down[[i]]
    
  )
  
}


# ----------------------------------------------------------
# 05.14 Add KEGG results to supplementary workbook
# ----------------------------------------------------------

addWorksheet(
  
  wb,
  
  "KEGG_Up"
  
)

writeData(
  
  wb,
  
  "KEGG_Up",
  
  KEGG_up[[1]]
  
)


addWorksheet(
  
  wb,
  
  "KEGG_Down"
  
)

writeData(
  
  wb,
  
  "KEGG_Down",
  
  KEGG_down[[1]]
  
)


# ----------------------------------------------------------
# 05.15 Add Reactome results to supplementary workbook
# ----------------------------------------------------------

addWorksheet(
  
  wb,
  
  "Reactome_Up"
  
)

writeData(
  
  wb,
  
  "Reactome_Up",
  
  Reactome_up[[1]]
  
)


addWorksheet(
  
  wb,
  
  "Reactome_Down"
  
)

writeData(
  
  wb,
  
  "Reactome_Down",
  
  Reactome_down[[1]]
  
)


# ----------------------------------------------------------
# 05.16 Save supplementary enrichment workbook
# ----------------------------------------------------------

saveWorkbook(
  
  wb,
  
  file.path(
    
    project_dir,
    
    "10_tables",
    
    "Supplementary_Enrichment.xlsx"
    
  ),
  
  overwrite = TRUE
  
)


# ----------------------------------------------------------
# 05.17 Save GO R objects
# ----------------------------------------------------------

saveRDS(
  
  GO_up,
  
  file.path(
    
    project_dir,
    
    "09_enrichment",
    
    "GO",
    
    "GO_Up.rds"
    
  )
  
)

saveRDS(
  
  GO_down,
  
  file.path(
    
    project_dir,
    
    "09_enrichment",
    
    "GO",
    
    "GO_Down.rds"
    
  )
  
)


# ----------------------------------------------------------
# 05.18 Save KEGG R objects
# ----------------------------------------------------------

saveRDS(
  
  KEGG_up,
  
  file.path(
    
    project_dir,
    
    "09_enrichment",
    
    "KEGG",
    
    "KEGG_Up.rds"
    
  )
  
)

saveRDS(
  
  KEGG_down,
  
  file.path(
    
    project_dir,
    
    "09_enrichment",
    
    "KEGG",
    
    "KEGG_Down.rds"
    
  )
  
)


# ----------------------------------------------------------
# 05.19 Save Reactome R objects
# ----------------------------------------------------------

saveRDS(
  
  Reactome_up,
  
  file.path(
    
    project_dir,
    
    "09_enrichment",
    
    "Reactome",
    
    "Reactome_Up.rds"
    
  )
  
)

saveRDS(
  
  Reactome_down,
  
  file.path(
    
    project_dir,
    
    "09_enrichment",
    
    "Reactome",
    
    "Reactome_Down.rds"
    
  )
  
)


# ----------------------------------------------------------
# 05.20 Enrichment analysis summary
# ----------------------------------------------------------

cat(
  
  "\nGO enrichment completed.\n",
  
  "KEGG enrichment completed.\n",
  
  "Reactome enrichment completed.\n"
  
)


# ----------------------------------------------------------
# 05.21 Load packages for enrichment visualization
# ----------------------------------------------------------

library(ggplot2)
library(dplyr)
library(stringr)
library(readr)


# ----------------------------------------------------------
# 05.22 Create enrichment figure directories
# ----------------------------------------------------------

dir.create(
  
  file.path(
    project_dir,
    "09_enrichment",
    "Figures"
  ),
  
  recursive = TRUE,
  
  showWarnings = FALSE
  
)

dir.create(
  
  file.path(
    project_dir,
    "09_enrichment",
    "Figures",
    "PDF"
  ),
  
  recursive = TRUE,
  
  showWarnings = FALSE
  
)

dir.create(
  
  file.path(
    project_dir,
    "09_enrichment",
    "Figures",
    "EPS"
  ),
  
  recursive = TRUE,
  
  showWarnings = FALSE
  
)


# ----------------------------------------------------------
# 05.23 General enrichment processing function
# ----------------------------------------------------------

prepare_enrichment <- function(
    df,
    group
){
  
  df %>%
    
    filter(
      P.value < 0.05
    ) %>%
    
    mutate(
      
      Count = lengths(
        
        strsplit(
          as.character(Genes),
          ";"
        )
        
      ),
      
      Score = -log10(
        P.value
      ),
      
      Term = gsub(
        
        "\\s*\\(GO:\\d+\\)|GO:\\d+\\s*\\|\\s*",
        
        "",
        
        Term
        
      ),
      
      Group = group
      
    ) %>%
    
    arrange(
      
      P.value,
      
      Adjusted.P.value,
      
      desc(Combined.Score)
      
    )
  
}


# ----------------------------------------------------------
# 05.24 Prepare GO enrichment data
# ----------------------------------------------------------

go_up_all <- bind_rows(
  
  prepare_enrichment(
    GO_up$GO_Biological_Process_2026,
    "BP"
  ),
  
  prepare_enrichment(
    GO_up$GO_Molecular_Function_2026,
    "MF"
  ),
  
  prepare_enrichment(
    GO_up$GO_Cellular_Component_2026,
    "CC"
  )
  
)


go_down_all <- bind_rows(
  
  prepare_enrichment(
    GO_down$GO_Biological_Process_2026,
    "BP"
  ),
  
  prepare_enrichment(
    GO_down$GO_Molecular_Function_2026,
    "MF"
  ),
  
  prepare_enrichment(
    GO_down$GO_Cellular_Component_2026,
    "CC"
  )
  
)


# ----------------------------------------------------------
# 05.25 Create GO supplementary tables
# ----------------------------------------------------------

GO_Up_Table <- go_up_all %>%
  
  arrange(
    P.value,
    Adjusted.P.value
  ) %>%
  
  dplyr::select(
    
    Term,
    Overlap,
    Count,
    P.value,
    Adjusted.P.value,
    Odds.Ratio,
    Combined.Score,
    Genes
    
  )


GO_Down_Table <- go_down_all %>%
  
  arrange(
    P.value,
    Adjusted.P.value
  ) %>%
  
  dplyr::select(
    
    Term,
    Overlap,
    Count,
    P.value,
    Adjusted.P.value,
    Odds.Ratio,
    Combined.Score,
    Genes
    
  )


# ----------------------------------------------------------
# 05.26 Select representative GO terms
# ----------------------------------------------------------

# Biological Process - Upregulated

keep_BP_up <- c(
  
  "Phagocytosis",
  
  "Inflammatory Response",
  
  "Complement Receptor Mediated Signaling Pathway",
  
  "Leukocyte Cell-Cell Adhesion",
  
  "Integrin-mediated Signaling Pathway"
  
)


GO_BP_Up <- go_up_all %>%
  
  filter(
    
    Group == "BP",
    
    Term %in% keep_BP_up
    
  )


# Molecular Function - Upregulated

keep_MF_up <- c(
  
  "Complement Receptor Activity",
  
  "Immunoglobulin Receptor Activity",
  
  "Protein Kinase Binding",
  
  "Toll-like Receptor Binding",
  
  "Superoxide-generating NADPH Oxidase Activator Activity"
  
)


GO_MF_Up <- go_up_all %>%
  
  filter(
    
    Group == "MF",
    
    Term %in% keep_MF_up
    
  )


# Cellular Component - Upregulated

keep_CC_up <- c(
  
  "Secretory Granule Membrane",
  
  "Tertiary Granule",
  
  "Cytoplasmic Vesicle Membrane",
  
  "Actin Filament",
  
  "Focal Adhesion"
  
)


GO_CC_Up <- go_up_all %>%
  
  filter(
    
    Group == "CC",
    
    Term %in% keep_CC_up
    
  )


# ----------------------------------------------------------
# 05.27 Rebuild selected GO Up dataset
# ----------------------------------------------------------

go_up <- bind_rows(
  
  GO_BP_Up,
  
  GO_MF_Up,
  
  GO_CC_Up
  
)


# ----------------------------------------------------------
# 05.28 Select representative GO Down terms
# ----------------------------------------------------------

# Biological Process - Downregulated

keep_BP_down <- c(
  
  "Cholesterol Biosynthetic Process",
  
  "Regulation of Lysosome Size",
  
  "Semaphorin-plexin Signaling Pathway",
  
  "Monoatomic Cation Transmembrane Transport",
  
  "Cellular Response to ATP"
  
)


GO_BP_Down <- go_down_all %>%
  
  filter(
    
    Group == "BP",
    
    Term %in% keep_BP_down
    
  )


# Molecular Function - Downregulated

keep_MF_down <- c(
  
  "Outward Rectifier Potassium Channel Activity",
  
  "Neuropeptide Y Receptor Activity",
  
  "Monoatomic Ion Channel Activity",
  
  "G Protein-Coupled Purinergic Nucleotide Receptor Activity",
  
  "Neuropilin Binding"
  
)


GO_MF_Down <- go_down_all %>%
  
  filter(
    
    Group == "MF",
    
    Term %in% keep_MF_down
    
  )


# Cellular Component - Downregulated

keep_CC_down <- c(
  
  "Endoplasmic Reticulum Lumen",
  
  "Late Endosome Membrane",
  
  "Cytoplasmic Side of Lysosomal Membrane",
  
  "Calcium Channel Complex",
  
  "Intermediate-density Lipoprotein Particle"
  
)


GO_CC_Down <- go_down_all %>%
  
  filter(
    
    Group == "CC",
    
    Term %in% keep_CC_down
    
  )


# ----------------------------------------------------------
# 05.29 Rebuild selected GO Down dataset
# ----------------------------------------------------------

go_down <- bind_rows(
  
  GO_BP_Down,
  
  GO_MF_Down,
  
  GO_CC_Down
  
)


# ----------------------------------------------------------
# 05.30 Define publication theme
# ----------------------------------------------------------

theme_pub <- theme_classic(
  
  base_family = "Arial"
  
) +
  
  theme(
    
    axis.text = element_text(
      colour = "black"
    ),
    
    axis.title = element_text(
      colour = "black"
    ),
    
    legend.text = element_text(
      colour = "black"
    ),
    
    legend.title = element_text(
      colour = "black"
    )
    
  )


# ----------------------------------------------------------
# 05.31 Define GO plotting function
# ----------------------------------------------------------

plot_GO <- function(
    df
){
  
  df <- df %>%
    
    arrange(
      
      P.value,
      
      Adjusted.P.value,
      
      desc(Combined.Score)
      
    ) %>%
    
    mutate(
      
      Term = stringr::str_wrap(
        
        Term,
        
        width = 40
        
      ),
      
      Term = factor(
        
        Term,
        
        levels = rev(
          unique(Term)
        )
        
      )
      
    )
  
  
  ggplot(
    
    df,
    
    aes(
      
      x = Score,
      
      y = Term,
      
      size = Count
      
    )
    
  ) +
    
    geom_point(
      
      aes(
        fill = Score
      ),
      
      shape = 21,
      
      colour = "black",
      
      stroke = 0.3
      
    ) +
    
    scale_fill_gradient(
      
      low = "#542788",
      
      high = "#F1A340",
      
      name = expression(
        -log[10](italic(P))
      )
      
    ) +
    
    scale_size(
      
      range = c(
        3.5,
        8
      ),
      
      name = "Gene Count"
      
    ) +
    
    scale_x_continuous(
      
      expand = expansion(
        
        mult = c(
          0.02,
          0.05
        )
        
      )
      
    ) +
    
    labs(
      
      x = expression(
        -log[10](italic(P))
      ),
      
      y = NULL
      
    ) +
    
    theme_pub +
    
    theme(
      
      axis.text.y = element_text(
        
        size = 10,
        
        colour = "black"
        
      ),
      
      axis.text.x = element_text(
        
        size = 10,
        
        colour = "black"
        
      ),
      
      legend.key.height = unit(
        0.8,
        "cm"
      ),
      
      panel.grid.major = element_blank(),
      
      panel.grid.minor = element_blank()
      
    )
  
}


# ----------------------------------------------------------
# 05.32 Draw GO figures
# ----------------------------------------------------------

Fig_GO_Up <-
  
  plot_GO(
    go_up
  )


Fig_GO_Down <-
  
  plot_GO(
    go_down
  )


# ----------------------------------------------------------
# 05.33 Save GO figures
# ----------------------------------------------------------

figure_dir <- file.path(
  
  project_dir,
  
  "09_enrichment",
  
  "Figures"
  
)


if (!dir.exists(
  figure_dir
)) {
  
  dir.create(
    
    figure_dir,
    
    recursive = TRUE
    
  )
  
}


ggsave(
  
  filename = file.path(
    
    figure_dir,
    
    "Fig4_GO_Up.tiff"
    
  ),
  
  plot = Fig_GO_Up,
  
  dpi = 600,
  
  width = 9,
  
  height = 6,
  
  compression = "lzw"
  
)


ggsave(
  
  filename = file.path(
    
    figure_dir,
    
    "Fig5_GO_Down.tiff"
    
  ),
  
  plot = Fig_GO_Down,
  
  dpi = 600,
  
  width = 9,
  
  height = 6,
  
  compression = "lzw"
  
)


# ----------------------------------------------------------
# 05.34 Prepare KEGG enrichment data
# ----------------------------------------------------------

kegg_up_all <- prepare_enrichment(
  
  KEGG_up[[1]],
  
  "Up"
  
)


kegg_down_all <- prepare_enrichment(
  
  KEGG_down[[1]],
  
  "Down"
  
)


# ----------------------------------------------------------
# 05.35 Create KEGG supplementary tables
# ----------------------------------------------------------

KEGG_Up_Table <-
  
  kegg_up_all %>%
  
  arrange(
    
    P.value,
    
    Adjusted.P.value,
    
    desc(Combined.Score)
    
  ) %>%
  
  dplyr::select(
    
    Term,
    Overlap,
    Count,
    P.value,
    Adjusted.P.value,
    Odds.Ratio,
    Combined.Score,
    Genes
    
  )


KEGG_Down_Table <-
  
  kegg_down_all %>%
  
  arrange(
    
    P.value,
    
    Adjusted.P.value,
    
    desc(Combined.Score)
    
  ) %>%
  
  dplyr::select(
    
    Term,
    Overlap,
    Count,
    P.value,
    Adjusted.P.value,
    Odds.Ratio,
    Combined.Score,
    Genes
    
  )


# ----------------------------------------------------------
# 05.36 Select representative KEGG pathways
# ----------------------------------------------------------

keep_KEGG_up <- c(
  
  "NEUTROPHIL EXTRACELLULAR TRAP FORMATION",
  
  "PHAGOSOME",
  
  "FC GAMMA R-MEDIATED PHAGOCYTOSIS",
  
  "CHEMOKINE SIGNALING PATHWAY",
  
  "LEUKOCYTE TRANSENDOTHELIAL MIGRATION",
  
  "LIPID AND ATHEROSCLEROSIS",
  
  "RAP1 SIGNALING PATHWAY"
  
)


kegg_up <-
  
  kegg_up_all %>%
  
  filter(
    
    Term %in% keep_KEGG_up
    
  )


keep_KEGG_down <- c(
  
  "TERPENOID BACKBONE BIOSYNTHESIS",
  
  "CHOLESTEROL METABOLISM",
  
  "REGULATION OF LIPOLYSIS IN ADIPOCYTES",
  
  "GNRH SECRETION",
  
  "PROTEIN DIGESTION AND ABSORPTION",
  
  "ALDOSTERONE SYNTHESIS AND SECRETION"
  
)


kegg_down <-
  
  kegg_down_all %>%
  
  filter(
    
    Term %in% keep_KEGG_down
    
  )


# ----------------------------------------------------------
# 05.37 Prepare Reactome enrichment data
# ----------------------------------------------------------

reactome_up_all <- prepare_enrichment(
  
  Reactome_up[[1]],
  
  "Up"
  
)


reactome_down_all <- prepare_enrichment(
  
  Reactome_down[[1]],
  
  "Down"
  
)


# ----------------------------------------------------------
# 05.38 Create Reactome supplementary tables
# ----------------------------------------------------------

Reactome_Up_Table <-
  
  reactome_up_all %>%
  
  arrange(
    
    P.value,
    
    Adjusted.P.value,
    
    desc(Combined.Score)
    
  ) %>%
  
  dplyr::select(
    
    Term,
    Overlap,
    Count,
    P.value,
    Adjusted.P.value,
    Odds.Ratio,
    Combined.Score,
    Genes
    
  )


Reactome_Down_Table <-
  
  reactome_down_all %>%
  
  arrange(
    
    P.value,
    
    Adjusted.P.value,
    
    desc(Combined.Score)
    
  ) %>%
  
  dplyr::select(
    
    Term,
    Overlap,
    Count,
    P.value,
    Adjusted.P.value,
    Odds.Ratio,
    Combined.Score,
    Genes
    
  )


# ----------------------------------------------------------
# 05.39 Select representative Reactome pathways
# ----------------------------------------------------------

keep_Reactome_up <- c(
  
  "Neutrophil Degranulation",
  
  "RHO GTPases Activate NADPH Oxidases",
  
  "ROS and RNS Production in Phagocytes",
  
  "Signaling by Interleukins",
  
  "Adaptive Immune System",
  
  "Antigen processing-Cross Presentation",
  
  "Cytokine Signaling in Immune System"
  
)


reactome_up <-
  
  reactome_up_all %>%
  
  filter(
    
    Term %in% keep_Reactome_up
    
  )


keep_Reactome_down <- c(
  
  "Cholesterol Biosynthesis",
  
  "Activation of Gene Expression by SREBF (SREBP)",
  
  "Plasma Lipoprotein Remodeling",
  
  "Collagen Formation",
  
  "NCAM Signaling for Neurite Out-Growth",
  
  "Tandem Pore Domain Potassium Channels",
  
  "Signaling by GPCR"
  
)


reactome_down <-
  
  reactome_down_all %>%
  
  filter(
    
    Term %in% keep_Reactome_down
    
  )


# ----------------------------------------------------------
# 05.40 Define KEGG plotting function
# ----------------------------------------------------------

plot_KEGG <- function(
    df
){
  
  df <- df %>%
    
    arrange(
      
      P.value,
      
      Adjusted.P.value,
      
      desc(Combined.Score)
      
    ) %>%
    
    mutate(
      
      Term = stringr::str_wrap(
        
        Term,
        
        width = 40
        
      ),
      
      Term = factor(
        
        Term,
        
        levels = rev(
          unique(Term)
        )
        
      )
      
    )
  
  
  ggplot(
    
    df,
    
    aes(
      
      x = Score,
      
      y = Term,
      
      size = Count
      
    )
    
  ) +
    
    geom_point(
      
      aes(
        fill = Score
      ),
      
      shape = 21,
      
      colour = "black",
      
      stroke = 0.3
      
    ) +
    
    scale_fill_gradient(
      
      low = "#542788",
      
      high = "#F1A340",
      
      name = expression(
        -log[10](italic(P))
      )
      
    ) +
    
    scale_size(
      
      range = c(
        3.5,
        8
      ),
      
      name = "Gene Count"
      
    ) +
    
    scale_x_continuous(
      
      expand = expansion(
        
        mult = c(
          0.02,
          0.05
        )
        
      )
      
    ) +
    
    labs(
      
      x = expression(
        -log[10](italic(P))
      ),
      
      y = NULL
      
    ) +
    
    theme_pub +
    
    theme(
      
      axis.text.y = element_text(
        
        size = 10,
        
        colour = "black"
        
      ),
      
      axis.text.x = element_text(
        
        size = 10,
        
        colour = "black"
        
      ),
      
      legend.key.height = unit(
        0.8,
        "cm"
      ),
      
      panel.grid.major = element_blank(),
      
      panel.grid.minor = element_blank()
      
    )
  
}


# ----------------------------------------------------------
# 05.41 Define Reactome plotting function
# ----------------------------------------------------------

plot_Reactome <- function(
    df
){
  
  df <- df %>%
    
    arrange(
      
      P.value,
      
      Adjusted.P.value,
      
      desc(Combined.Score)
      
    ) %>%
    
    mutate(
      
      Score = pmin(
        Score,
        20
      ),
      
      Term = factor(
        
        Term,
        
        levels = rev(
          unique(Term)
        )
        
      )
      
    )
  
  
  ggplot(
    
    df,
    
    aes(
      
      x = Score,
      
      y = Term,
      
      size = Count
      
    )
    
  ) +
    
    geom_point(
      
      aes(
        fill = Score
      ),
      
      shape = 21,
      
      colour = "black",
      
      stroke = 0.3
      
    ) +
    
    scale_fill_gradient(
      
      low = "#1B9E77",
      
      high = "#D95F02",
      
      name = expression(
        -log[10](italic(P))
      )
      
    ) +
    
    scale_size(
      
      range = c(
        3.5,
        8
      ),
      
      name = "Gene Count"
      
    ) +
    
    scale_x_continuous(
      
      expand = expansion(
        
        mult = c(
          0.02,
          0.05
        )
        
      )
      
    ) +
    
    labs(
      
      x = expression(
        -log[10](italic(P))
      ),
      
      y = NULL
      
    ) +
    
    theme_pub +
    
    theme(
      
      axis.text.y = element_text(
        
        size = 10,
        
        colour = "black"
        
      ),
      
      axis.text.x = element_text(
        
        size = 10,
        
        colour = "black"
        
      ),
      
      legend.key.height = unit(
        0.8,
        "cm"
      ),
      
      panel.grid.major = element_blank(),
      
      panel.grid.minor = element_blank()
      
    )
  
}


# ----------------------------------------------------------
# 05.42 Generate KEGG figures
# ----------------------------------------------------------

Fig_KEGG_Up <-
  
  plot_KEGG(
    kegg_up
  )


Fig_KEGG_Down <-
  
  plot_KEGG(
    kegg_down
  )


# ----------------------------------------------------------
# 05.43 Generate Reactome figures
# ----------------------------------------------------------

Fig_Reactome_Up <-
  
  plot_Reactome(
    reactome_up
  )


Fig_Reactome_Down <-
  
  plot_Reactome(
    reactome_down
  )


# ----------------------------------------------------------
# 05.44 Define figure and supplementary directories
# ----------------------------------------------------------

figure_dir <- file.path(
  
  project_dir,
  
  "09_enrichment",
  
  "Figures"
  
)


table_dir <- file.path(
  
  project_dir,
  
  "09_enrichment",
  
  "Supplementary"
  
)


if (!dir.exists(
  figure_dir
)) {
  
  dir.create(
    
    figure_dir,
    
    recursive = TRUE
    
  )
  
}


if (!dir.exists(
  table_dir
)) {
  
  dir.create(
    
    table_dir,
    
    recursive = TRUE
    
  )
  
}


# ----------------------------------------------------------
# 05.45 Save KEGG figures
# ----------------------------------------------------------

ggsave(
  
  file.path(
    
    figure_dir,
    
    "Fig6_KEGG_Up.tiff"
    
  ),
  
  Fig_KEGG_Up,
  
  dpi = 600,
  
  width = 9,
  
  height = 5,
  
  compression = "lzw"
  
)


ggsave(
  
  file.path(
    
    figure_dir,
    
    "Fig7_KEGG_Down.tiff"
    
  ),
  
  Fig_KEGG_Down,
  
  dpi = 600,
  
  width = 9,
  
  height = 5,
  
  compression = "lzw"
  
)


# ----------------------------------------------------------
# 05.46 Save Reactome figures
# ----------------------------------------------------------

ggsave(
  
  file.path(
    
    figure_dir,
    
    "Fig8_Reactome_Up.tiff"
    
  ),
  
  Fig_Reactome_Up,
  
  dpi = 600,
  
  width = 9,
  
  height = 5,
  
  compression = "lzw"
  
)


ggsave(
  
  file.path(
    
    figure_dir,
    
    "Fig9_Reactome_Down.tiff"
    
  ),
  
  Fig_Reactome_Down,
  
  dpi = 600,
  
  width = 9,
  
  height = 5,
  
  compression = "lzw"
  
)


# ----------------------------------------------------------
# 05.47 Export supplementary enrichment tables
# ----------------------------------------------------------

library(
  writexl
)


write_xlsx(
  
  GO_Up_Table,
  
  file.path(
    
    table_dir,
    
    "Supplementary_Table_S5_GO_Up.xlsx"
    
  )
  
)


write_xlsx(
  
  GO_Down_Table,
  
  file.path(
    
    table_dir,
    
    "Supplementary_Table_S6_GO_Down.xlsx"
    
  )
  
)


write_xlsx(
  
  KEGG_Up_Table,
  
  file.path(
    
    table_dir,
    
    "Supplementary_Table_S7_KEGG_Up.xlsx"
    
  )
  
)


write_xlsx(
  
  KEGG_Down_Table,
  
  file.path(
    
    table_dir,
    
    "Supplementary_Table_S8_KEGG_Down.xlsx"
    
  )
  
)


write_xlsx(
  
  Reactome_Up_Table,
  
  file.path(
    
    table_dir,
    
    "Supplementary_Table_S9_Reactome_Up.xlsx"
    
  )
  
)


write_xlsx(
  
  Reactome_Down_Table,
  
  file.path(
    
    table_dir,
    
    "Supplementary_Table_S10_Reactome_Down.xlsx"
    
  )
  
)
