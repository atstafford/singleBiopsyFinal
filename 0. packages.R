#devtools::install_github("psyteachr/introdataviz")

Packages <- c(
              
              "Hmisc", # distribution of subclonal
              "devtools" #for paired violin
              ) 

Packages <- c("tidyverse","readxl","cowplot","ggrepel","ggpubr", # generic
              "reshape2", # for melting in HM
              "DescTools", # for overlap with telo/centro
              "dplyr", # for kw test
              "GenomicRanges", "rtracklayer", # for hg conversions
              "msigdbr", # for hallmark gene lists
              "SciViews", "vegan", # for shannon calculations
              "pvclust", # clustering
              "fastDummies","betareg","car","frmselection", # building model
              "caret", #LOOCV
              "BSgenome", # for hg19 coordinate in chromomap
              "chromoMap", # visualising bin location
              "org.Hs.eg.db", # swtiching between gene IDS
              "limma","rrvgo", # GO
              "survival","survminer", # for survival
              "TCGAbiolinks", # pull tcga rna data
              "TCGAutils", # for TCGA UUID barcode translator
              "DESeq2","powerjoin" # for RNA expression
              )

lapply(Packages, library, character.only = TRUE)


