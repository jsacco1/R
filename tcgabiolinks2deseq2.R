rm(list = ls()) #clear all environmental variables

# Title: TCGAbiolinks and DESeq2 Analysis ####

# install with BiocManager::install(c("TCGAbiolinks", "DESeq2"))
# I personally recommend not installing and running this analysis via RStudio 
# conda environment.

library(TCGAbiolinks)
library(DESeq2)

# Query GDC ####
# internet connection required
# pick a TCGA project. for example: "TCGA-BRCA", "TCGA-GBM", "TCGA-OV".
# see package docs for complete list of valid projects.

proj <- "TCGA-BRCA"
query <- GDCquery(
  project = proj,
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts"
)

# Download Data ####

GDCdownload(query)
data <- GDCprepare(query)

# Data Processing ####

# remove NA values
data <- data[,!is.na(data$paper_IDH.status)]

ddsSE <- DESeqDataSet(data, design = ~ paper_IDH.status)

# Data Processing ####

# rename first element in assays to "counts"
# coerce counts to integer datatype

keep <- rowSums(counts(ddsSE)) >= 10
ddsSE <- ddsSE[keep,]
ddsSE <- DESeq(ddsSE)

resultsNames(ddsSE)

# check levels
ddsSE$paper_IDH.status
# Levels: Mutant WT

# Summarized Results ####

res <- results(ddsSE, name = "paper_IDH.status_WT_vs_Mutant")
res_df <- as.data.frame(res)
summary(res)

# out of 46017 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 3331, 7.2%
# LFC < 0 (down)     : 2739, 6%
# outliers [1]       : 0, 0%
# low counts [2]     : 12492, 27%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Errata ####

# Reviewing the TCGAbiolinks Github issues page, I found that 
# installing the following packages may help if the following 
# error arises:

# Error in value[[3L]](cond) : 
# GDC server down, try to use this package later

# Packages:
# install.packages("GlobalOptions")
# BiocManager::install('Bioconductor/GenomicDataCommons')
# BiocManager::install('BioinformaticsFMRP/TCGAbiolinks',ref = "GenomicDataCommons")
# devtools::install_github(repo = "grimbough/biomaRt",dependencies = T)
# devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")

# Note that only some of the above may successfully install.

# References ####

# https://rpubs.com/tiagochst/TCGAbiolinks_to_DESEq2

# Session Info ####
devtools::session_info()
