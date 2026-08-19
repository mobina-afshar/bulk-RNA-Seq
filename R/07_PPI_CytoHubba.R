############################################################
# 07. Protein-Protein Interaction (PPI) Network Analysis
############################################################

# This step was performed using STRING and Cytoscape.
#
# Differentially expressed genes were used to construct
# protein-protein interaction networks, followed by
# network visualization and hub-gene ranking using
# Cytoscape/CytoHubba.
#
# The resulting CytoHubba rankings were imported in the
# subsequent consensus hub-gene analysis.


# ----------------------------------------------------------
# 07.1 PPI network construction
# ----------------------------------------------------------

# PPI networks were constructed using STRING and
# visualized in Cytoscape.
#
# The resulting network files are not included in this
# public repository because they were generated externally
# and are not part of the R analysis script.


# ----------------------------------------------------------
# 07.2 Hub-gene identification
# ----------------------------------------------------------

# Hub genes were ranked using four CytoHubba algorithms:
#
#   - MCC
#   - Degree
#   - MNC
#   - EPC
#
# The resulting ranking tables were exported from
# Cytoscape and used as input for the consensus
# hub-gene analysis.


# ----------------------------------------------------------
# 07.3 Downstream analysis
# ----------------------------------------------------------

# The exported CytoHubba results are imported and
# integrated in:
#
#   08_consensus_hub_genes.R
#
# Consensus hub genes are defined based on their
# recurrence across the four ranking algorithms.
