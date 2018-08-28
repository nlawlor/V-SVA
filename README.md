# Iteratively Adjusted Surrogate Variable Analysis (*IA-SVA*)

### A shiny application for interactive analysis and visualization of transcriptomic data

Authors: *Donghyung Lee*, *Anthony Cheng*, *Nathan Lawlor*, *Mohan Bolisetty*, and *Duygu Ucar*

App Maintainer: *Nathan Lawlor*

***

### Introduction

*IA-SVA* is a statistical framework to uncover hidden sources of variation even when these sources are correlated with the biological variable of interest. *IA-SVA* provides a flexible methodology to i) identify a hidden factor for unwanted heterogeneity while adjusting for all known factors; ii) test the significance of the putative hidden factor for explaining the variation in the data; and iii), if significant, use the estimated factor as an additional known factor in the next iteration to uncover further hidden factors.

This shiny app reads a matrix of gene expression data (rows as genes, columns as samples) and a matrix of sample metadata (rows as samples, columns as different traits) and provides a suite of methods for feature extraction, gene set enrichment, and visualization of transcriptomic data.

***

### Tutorial/Guide

For an in depth guide of how to use this shiny app, please see the Word Document at: https://github.com/nlawlor/iasva_shiny/blob/master/Data/IASVA_Shiny_App_Manual.docx

### Required Input Files 

1. The app requires first a matrix of gene expression data (rows as genes, columns as samples). This file may be formatted as a tab-delimited text file, .CSV file, or .Rds object.  

An example 10X Genomics single cell RNA-seq expression matrix obtained from [*Kang et al. 2017*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5784859/) which contains a subset of 200 human peripheral blood mononuclear cells (PBMCs) 
is provided at: https://github.com/nlawlor/iasva_shiny/blob/master/Data/test.exp.txt and can be used to test the app.

This matrix, when loaded into R:

```R, echo=FALSE, message=FALSE, eval=TRUE
df <- read.delim("Data/test.exp.txt", header = T, check.names = F, stringsAsFactors = F, row.names = 1)
df[1:4, 1:4]
```

Looks like this: ![](https://github.com/nlawlor/iasva_shiny/blob/master/img/exp.matrix.png)

2. The app secondly requires a matrix of sample metadata (rows as sample names, columns as traits). Like the expression matrix, this file may be formatted as a tab-delimited text file, .CSV file, or .Rds object.  

An example metadata file for the above expression matrix can be found at: https://github.com/nlawlor/iasva_shiny/blob/master/Data/test.metadata.txt.

This matrix, when loaded into R:

```R, echo=FALSE, message=FALSE, eval=TRUE
meta <- read.delim("Data/test.metadata.txt", header = T, check.names = F, stringsAsFactors = F, row.names = 1)
meta[1:4, ]
```
Looks like this: ![](https://github.com/nlawlor/iasva_shiny/blob/master/img/metadata.png)

**Please Note: The sample identifiers contained in the rows of the metadata matrix should be the same as the identifiers provided in the columns of the expression matrix**

***

### Overview of features included in IA-SVA shiny

* **Identify hidden sources of heterogeneity in transcriptomic data:**

![](https://github.com/nlawlor/iasva_shiny/blob/master/img/sv.plots.png)

* **Discover marker genes associated with a variable of interest:**

![](https://github.com/nlawlor/iasva_shiny/blob/master/img/marker.genes.png)

* **Determine molecular pathways associated with the identified marker genes:**

![](https://github.com/nlawlor/iasva_shiny/blob/master/img/pathway.analysis.png)

* **Clustering and interactive visualization of data:**

![](https://github.com/nlawlor/iasva_shiny/blob/master/img/tsne.gif)

***

### Availability

* The app is currently hosted on **shinyapps.io** and available for use [here](https://nlawlor.shinyapps.io/IASVA_Shiny_08_13_2018/).

* **IA-SVA** is also available as an R-package that may be downloaded from **Bioconductor** [here](https://www.bioconductor.org/packages/devel/bioc/html/iasva.html), and the source code is available [here](https://github.com/UcarLab/iasva)

### Other Resources

* The preprint manuscript describing this software can be found on **bioRxiv** [here](https://www.biorxiv.org/content/early/2018/04/24/151217)
