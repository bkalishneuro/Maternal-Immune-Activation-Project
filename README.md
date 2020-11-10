# Maternal-Immune-Activation-Project
This repository organizes custom code used by Brian Kalish and Benjamin Finander to analyze data for the manuscript, "Maternal Immune Activation in Mice Disrupts Proteostasis in the Fetal Brain"

The code is organized by analysis: Ribosome Profiling and scRNA-seq

Notes:
-In instances where the same bioinformatic analysis was applied several times (i.e. multiple cell clusters), one example script will be posted

-For the Ribosome Profiling Scripts, Ribo-seq and bulk RNA-seq are processed separately to generate fastq counts files from fastq files, and then processed in the same manner for differential gene expression

-For the scRNA-seq scripts, the order of operations was as follows:


  Bcl files were converted to fastq files (general script)
  
  Fastq files were mapped using standard indrops procedure (specific yaml files and batch scripts for E14 and E18)
  
  Seurat objects were generated (specific scripts for E14 and E18)
  
  Undesired doublet cells and clusters were filtered (general script)
  
  After reclustering, seurat objects were divided into groups and subclustered (specific scripts for all)
  
  Normalized Counts objects were generated for all new seurat objects (general script)
  
  Differential gene expression was performed in EdgeR (general script)
