---
#title: 'RNA-Seq with knockdown'
#author: "James Sacco"
#date: "`r Sys.Date()`"
#output: 

# clear
rm(list=ls())

# load modules
library(GEOquery)
library(exprso)
studyID <- 'GSE159049'
gse = getGEO(GEO = studyID)
gse[[1]]

# get dependencies
sessionInfo()
