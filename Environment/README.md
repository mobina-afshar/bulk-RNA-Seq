# RNA-seq Analysis of Granulosa Cells in PCOS

This repository contains the computational analysis workflow used to investigate transcriptomic alterations in granulosa cells from patients with polycystic ovary syndrome (PCOS).

The workflow includes RNA-seq preprocessing, read alignment, gene-level quantification, differential expression analysis, functional enrichment, gene set enrichment analysis, protein-protein interaction analysis, consensus hub-gene identification, and lncRNA–hub gene correlation analysis.

## Project Overview

Polycystic ovary syndrome (PCOS) is a complex reproductive and metabolic disorder. This project investigates transcriptional alterations in granulosa cells associated with PCOS using publicly available transcriptomic data.

The analysis workflow was implemented using Linux shell scripts and R.

## Workflow

The analysis consists of the following major steps:

1. RNA-seq preprocessing and quality control
2. Reference genome alignment
3. Gene-level read quantification
4. Differential expression analysis
5. Differentially expressed gene annotation
6. Functional enrichment analysis
7. Gene Set Enrichment Analysis (GSEA)
8. Protein-protein interaction (PPI) network analysis
9. Consensus hub-gene identification
10. lncRNA–hub gene correlation analysis

## Repository Structure

```text
.
├── README.md
│
├── scripts/
│   ├── 01_preprocessing.sh
│   ├── 02_alignment.sh
│   └── 03_featureCounts.sh
│
├── R/
│   ├── 01_DESeq2.R
│   ├── 02_annotation_DEGs.R
│   ├── 03_visualization.R
│   ├── 04_lncRNA_analysis.R
│   ├── 05_functional_enrichment.R
│   ├── 06_GSEA.R
│   ├── 07_PPI_CytoHubba.R
│   ├── 08_consensus_hub_genes.R
│   └── 09_lncRNA_hub_correlation.R
│
└── environment/
    ├── README.md
    └── rnaseq_env.yml
``` 
## Computational Workflow
# Environment

The `rnaseq_env.yml` file contains the Conda environment specification used for the RNA-seq preprocessing and analysis workflow.

The environment can be recreated using:

```bash
conda env create -f rnaseq_env.yml
```
After creation, activate the environment with:
```bash
conda activate rnaseq
```
### 1. RNA-seq preprocessing

Raw sequencing reads were processed using Linux shell scripts.

The preprocessing workflow includes quality control and read preprocessing prior to alignment.

### 2. Read alignment

Processed reads were aligned to the human reference genome using STAR.

Alignment files were processed using SAMtools, and relevant alignment quality metrics were evaluated.

### 3. Gene-level quantification

Aligned reads were assigned to genes using featureCounts to generate a gene-level count matrix for downstream differential expression analysis.

### 4. Differential expression analysis

Differential expression analysis was performed in R using DESeq2.

The analysis included:

- count filtering
- normalization
- variance-stabilizing transformation (VST)
- principal component analysis (PCA)
- differential expression testing
- annotation of differentially expressed genes

### 5. Functional enrichment analysis

Functional enrichment analysis was performed using Enrichr.

The analysis included:

- Gene Ontology Biological Process (GO-BP)
- Gene Ontology Molecular Function (GO-MF)
- Gene Ontology Cellular Component (GO-CC)
- KEGG pathway enrichment
- Reactome pathway enrichment

Enrichment results were exported as tables and visualized using R.

### 6. Gene Set Enrichment Analysis

GSEA was performed using ranked gene-level statistics.

Three gene-set collections were evaluated:

- KEGG
- Reactome
- MSigDB Hallmark gene sets

GSEA results were exported and visualized using `clusterProfiler`, `ReactomePA`, `msigdbr`, and `enrichplot`.

### 7. PPI network analysis

Protein-protein interaction networks were constructed using STRING and visualized in Cytoscape.

Hub genes were ranked using four CytoHubba algorithms:

- Maximal Clique Centrality (MCC)
- Degree
- Maximum Neighborhood Component (MNC)
- Edge Percolated Component (EPC)

### 8. Consensus hub-gene identification

The rankings obtained from the four CytoHubba algorithms were integrated to identify consensus hub genes.

Genes detected by at least three of the four ranking algorithms were considered consensus hub genes.

An UpSet plot was generated to visualize the overlap among ranking methods.

### 9. lncRNA–hub gene correlation analysis

Expression relationships between differentially expressed lncRNAs and consensus hub genes were evaluated using Spearman correlation.

Correlation significance was assessed using P values, followed by Benjamini–Hochberg false discovery rate (FDR) correction.

Strong correlations were further examined using:

- absolute Spearman correlation coefficient ≥ 0.80
- P value < 0.05

Significant lncRNA–hub gene relationships were exported for downstream network visualization in Cytoscape.

## Software and R Packages

The analysis uses a combination of command-line bioinformatics tools and R packages.

### Major command-line and network-analysis tools

- FastQC
- Trim Galore
- STAR
- SAMtools
- featureCounts
- STRING
- Cytoscape
- CytoHubba

### Major R packages

- DESeq2
- clusterProfiler
- ReactomePA
- msigdbr
- enrichR
- enrichplot
- org.Hs.eg.db
- Hmisc
- ComplexHeatmap
- circlize
- ggplot2
- dplyr
- readr
- writexl
- openxlsx
- UpSetR

## Reproducibility

The repository provides the analysis scripts used for the computational workflow.

Raw sequencing data and intermediate analysis results are not included in this repository.

Dataset identifiers and complete methodological details should be consulted together with the associated manuscript.

## Data Availability

The analysis was performed using publicly available transcriptomic data.

Raw data are not redistributed in this repository.


## Contact

For questions regarding the analysis workflow or code, please open an issue in this repository.
