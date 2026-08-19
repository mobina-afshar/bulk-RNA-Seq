############################################################
# 08. Consensus Hub Gene Identification
############################################################


# ----------------------------------------------------------
# 08.1 Load required packages
# ----------------------------------------------------------

library(
  dplyr
)

library(
  readr
)

library(
  UpSetR
)

library(
  ggplot2
)


# ----------------------------------------------------------
# 08.2 Import CytoHubba results
# ----------------------------------------------------------

MCC <- read_csv(
  
  file.path(
    
    project_dir,
    
    "12_network_analysis",
    
    "string_interactions_short.tsv_MCC_top20_and_expanded default node.csv"
    
  )
  
)


Degree <- read_csv(
  
  file.path(
    
    project_dir,
    
    "12_network_analysis",
    
    "string_interactions_short.tsv_Degree_top20_and_expanded default node.csv"
    
  )
  
)


MNC <- read_csv(
  
  file.path(
    
    project_dir,
    
    "12_network_analysis",
    
    "string_interactions_short.tsv_MNC_top20_and_expanded default node.csv"
    
  )
  
)


EPC <- read_csv(
  
  file.path(
    
    project_dir,
    
    "12_network_analysis",
    
    "string_interactions_short.tsv_EPC_top20_and_expanded default node.csv"
    
  )
  
)


# ----------------------------------------------------------
# 08.3 Extract hub genes from each ranking method
# ----------------------------------------------------------

hub_lists <- list(
  
  MCC = MCC$name,
  
  Degree = Degree$name,
  
  MNC = MNC$name,
  
  EPC = EPC$name
  
)


# ----------------------------------------------------------
# 08.4 Calculate hub-gene occurrence
# ----------------------------------------------------------

hub_occurrence <- table(
  
  unlist(
    hub_lists
  )
  
)

hub_occurrence <- data.frame(
  
  Gene_Symbol = names(
    hub_occurrence
  ),
  
  Occurrence = as.numeric(
    hub_occurrence
  )
  
)


# ----------------------------------------------------------
# 08.5 Sort hub genes by occurrence
# ----------------------------------------------------------

hub_occurrence <- hub_occurrence %>%
  
  arrange(
    
    desc(
      Occurrence
    ),
    
    Gene_Symbol
    
  )


# ----------------------------------------------------------
# 08.6 Identify consensus hub genes
# ----------------------------------------------------------

# Consensus hub genes are defined as genes detected
# by at least three of the four CytoHubba algorithms.

Consensus_Hubs <- hub_occurrence %>%
  
  filter(
    
    Occurrence >= 3
    
  )


# ----------------------------------------------------------
# 08.7 Display consensus hub genes
# ----------------------------------------------------------

print(
  Consensus_Hubs
)

cat(
  
  "\nConsensus Hub Genes:",
  
  nrow(
    Consensus_Hubs
  ),
  
  "\n"
  
)


# ----------------------------------------------------------
# 08.8 Save hub-gene occurrence results
# ----------------------------------------------------------

write_csv(
  
  hub_occurrence,
  
  file.path(
    
    project_dir,
    
    "12_network_analysis",
    
    "Hub_Gene_Occurrence.csv"
    
  )
  
)


# ----------------------------------------------------------
# 08.9 Save consensus hub genes
# ----------------------------------------------------------

write_csv(
  
  Consensus_Hubs,
  
  file.path(
    
    project_dir,
    
    "12_network_analysis",
    
    "Consensus_Hub_Genes.csv"
    
  )
  
)


# ----------------------------------------------------------
# 08.10 Prepare UpSet plot data
# ----------------------------------------------------------

hub_sets <- list(
  
  MCC = MCC$name,
  
  Degree = Degree$name,
  
  MNC = MNC$name,
  
  EPC = EPC$name
  
)


upset_data <- fromList(
  hub_sets
)


# ----------------------------------------------------------
# 08.11 Generate UpSet plot
# ----------------------------------------------------------

tiff(
  
  filename = file.path(
    
    project_dir,
    
    "12_network_analysis",
    
    "Figures",
    
    "Consensus_Hub_UpSet.tiff"
    
  ),
  
  width = 7,
  
  height = 5,
  
  units = "in",
  
  res = 600,
  
  compression = "lzw"
  
)


upset(
  
  upset_data,
  
  sets = c(
    
    "MCC",
    "Degree",
    "MNC",
    "EPC"
    
  ),
  
  order.by = "freq",
  
  sets.bar.color = "grey40",
  
  main.bar.color = "steelblue",
  
  text.scale = 1.3
  
)


dev.off()


# ----------------------------------------------------------
# 08.12 Prepare consensus hub-gene expression information
# ----------------------------------------------------------

hub_expression <- Consensus_Hubs %>%
  
  left_join(
    
    results_annotation %>%
      
      dplyr::select(
        
        Gene_Symbol,
        
        log2FoldChange
        
      ),
    
    by = "Gene_Symbol"
    
  ) %>%
  
  mutate(
    
    Regulation = case_when(
      
      log2FoldChange > 0 ~ "Upregulated",
      
      log2FoldChange < 0 ~ "Downregulated",
      
      TRUE ~ "Unknown"
      
    )
    
  )


# ----------------------------------------------------------
# 08.13 Order consensus hub genes
# ----------------------------------------------------------

hub_expression <- hub_expression %>%
  
  arrange(
    
    Occurrence,
    
    abs(
      log2FoldChange
    )
    
  ) %>%
  
  mutate(
    
    Gene_Symbol = factor(
      
      Gene_Symbol,
      
      levels = Gene_Symbol
      
    )
    
  )


# ----------------------------------------------------------
# 08.14 Generate consensus hub-gene barplot
# ----------------------------------------------------------

Fig_Hub_Occurrence <-
  
  ggplot(
    
    hub_expression,
    
    aes(
      
      x = Occurrence,
      
      y = Gene_Symbol,
      
      fill = Regulation
      
    )
    
  ) +
  
  geom_col(
    
    width = 0.7
    
  ) +
  
  geom_text(
    
    aes(
      
      label = Occurrence
      
    ),
    
    hjust = -0.25,
    
    size = 4
    
  ) +
  
  scale_fill_manual(
    
    values = c(
      
      "Upregulated" = "#B40426",
      
      "Downregulated" = "#3B4CC0",
      
      "Unknown" = "grey60"
      
    ),
    
    name = NULL
    
  ) +
  
  scale_x_continuous(
    
    limits = c(
      
      0,
      
      4.5
      
    ),
    
    breaks = c(
      
      1,
      2,
      3,
      4
      
    )
    
  ) +
  
  labs(
    
    x = "Number of CytoHubba algorithms",
    
    y = NULL
    
  ) +
  
  theme_classic(
    
    base_family = "Arial",
    
    base_size = 12
    
  ) +
  
  theme(
    
    axis.text = element_text(
      
      colour = "black"
      
    ),
    
    axis.title.x = element_text(
      
      face = "bold"
      
    ),
    
    legend.position = "top"
    
  )


# ----------------------------------------------------------
# 08.15 Display consensus hub-gene barplot
# ----------------------------------------------------------

Fig_Hub_Occurrence


# ----------------------------------------------------------
# 08.16 Save consensus hub-gene barplot
# ----------------------------------------------------------

ggsave(
  
  filename = file.path(
    
    project_dir,
    
    "12_network_analysis",
    
    "Figures",
    
    "Consensus_Hub_Occurrence_Barplot.tiff"
    
  ),
  
  plot = Fig_Hub_Occurrence,
  
  width = 6,
  
  height = 7,
  
  dpi = 600,
  
  compression = "lzw"
  
)
