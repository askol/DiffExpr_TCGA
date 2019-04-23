Sex Differential Expression Analysis on The Cancer Genome Atlas Studies
==========

This repository contains scripts to perform differential expression analysis (DE) and evaluate and summarize the results within and across cancers.

Dependencies
------------

- All code is written in R.
- R libraries required include:
  - edgeR
  - limma
  - dplyr
  - survival
  - survminer
  - WGCNA
  - pscl
  - reshape2

- Code is designed to run on an HPC running PBS. 

Description of files
--------------------

filename                          |  description
----------------------------------|------------------------------------------------------------------------------------
README.md                         |  Text file (markdown format) description of the project and files.


R scripts:

filename                          |  description
----------------------------------|------------------------------------------------------------------------------------
DE_limma.r                        |  Runs sex DE. Takes TCGA study abreviation as argument. Write results to hard coded directory.
Run_DE_limma.r			  |  Loops through all TCGA studies, writes PBS script to run DE using DE_limma.r. Submits jobs to cluster.

Compare_DE_results.r		  |  Reads in results from DE analysis performed via Run_DE_limma.r. Creates files for and runs Gene Set Enrichment Analysis. Creates a file for the summarized output from the GSEA analysis.r
Compare_DE_results_funcs.r	  |  Contains the functions used to run GSEA, collect, and summarize the data.
Run_GSEA_analysis.r		  |  Loops through all TCGA studies, writes PBS script to run Compare_DE_results.r for each TCGA study. Submits jobs to cluster.

Summary.r			  |  Performs a number of analyses to examine the DE results, includes plots to summarize the number of DE genes per study, and the distribution of the number of DE genes per chromosome. Plot to compare, across studies, the distribution of DE genes across the X chromosome. Identified DE genes that are shared across multiple cancers, and between cancers and GTEx tissues. Etc.
Summary_funcs.r			  |  Functions to perform the analyses and plotting described requested from Summary.r
