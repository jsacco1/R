# R bioinformatics


This repository contains R scripts to calculate and graph incidence rates and allele frequencies. 

I have also included a run-through of a limma tutorial. Limma is a BioConductor package commonly used in RNA-Seq pipelines for differential expression analysis. 

### Query GDC

    Datasets can be accessed via the GDC repository at the following link:

https://portal.gdc.cancer.gov/repository?facetTab=files

### Query GDC for RNA-Seq data

    Filters:
    Data Category: transcriptome profiling
    Data Type: Gene Expression Quantification
    Experimental Strategy: RNA-Seq
    Workflow Type: HTSeq - Counts OR STAR - Counts
    Access: open

DGE analysis with DESeq2, edgeR, etc. requires raw counts, not FPKM. 

### Query GDC for WXS or Targeted Sequencing

    Filters:
    Data Type: Masked Somatic Mutation
    Experimental Strategy: WXS OR Targeted Sequencing
    Access: open
    
### RNA-Seq Data Formats
    raw counts (txt)
    RSEM (available from TCGA RNA-Seq V2; can be used with edgeR)

### Helpful links:
    Question: Interpreting TCGA .rsem.genes.results and .rsem.genes.normalized_results files.
    https://www.biostars.org/p/106127/ 
    
    Importing transcript abundance with tximport: RSEM
    https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#RSEM
    
