# PGG_DB
Pseudogene-gene network database and similarity matrix (beta version)

This directory contains the Supplementary Materials of our PSB 2018 paper.

Citation: Johnson TS, Li S, Kho JR, Huang K, Zhang Y. Network analysis of pseudogene-gene relationships: from pseudogene evolution to their functional potentials. Pac Symp Biocomput. 2018;23:536-547.

Access to paper: [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/29218912) [World Scientific](http://www.worldscientific.com/doi/abs/10.1142/9789813235533_0049) [PSB2018](https://psb.stanford.edu/psb-online/proceedings/psb18/)

# alignment_matrix.csv.zip
This file is an alignment matrix between the pseudogenes and gene families. Each row (except the last) is a gene family consensus sequence and each column is a pseudogene. The individual cells are the alignment score from ClustalW between the pseudogene sequence (column) and the consensus sequence (row). The last row is not a consensus sequence but the index (row number) of the highest scoring alignment in that column. This row number is the family to which the pseudogene was assigned. These new PGG families that are inclusive of pseudogenes and genes were then aligned within themselves (see below).

# pgAmats
This folder contains the complete adjacency matrix of each PGG family. Each family has its own csv file containing the sequence alignment score (of every pair of members in the family) from the CUDA-align algorithm where each row is a pseudogene/gene and each column is a pseudogene/gene. These were used to generate the MSTs used in our manuscript. There is also much more information contained within these files that was not fully reported in the manuscript due to page limit.

The zipped file containing all the adjacency matrices is also downloadable as **pgAmats.zip**.
