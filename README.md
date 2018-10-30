# Visual Iteratively Adjusted Surrogate Variable Analysis (*VIA-SVA*)

### An R Shiny application for detecting and annotating sources of variation in single cell RNA-seq data

Authors: *Donghyung Lee*, *Anthony Cheng*, *Nathan Lawlor*, *Mohan Bolisetty*, and *Duygu Ucar*

App Maintainer: *Nathan Lawlor*

App Location: https://nlawlor.shinyapps.io/VIASVA/

***

### Introduction

Single cell RNA-sequencing (scRNA-seq) is commonly used to define gene expression programs from individual cells. scRNA-seq expression data is highly spare and subject to substantial sources of experimental noise (“unwanted variation”) which can inhibit identification and analysis of unknown biological factors (“wanted variation”). Surrogate variable analysis (SVA) algorithms, such as *Iteratively Adjusted-Surrogate Variable Analysis (IA-SVA)*, have been developed to effectively estimate sources of hidden variation in expression data, however, annotating each factor and interpreting their underlying biology remains challenging. 

To facilitate the interpretation of detected hidden factors, we developed **Visual Iteratively Adjusted Surrogate Variable Analysis (VIA-SVA)**, an R Shiny application that provides a web-browser interface to *IA-SVA* for identification and annotation of hidden sources of variation in scRNA-seq data. This interactive framework includes tools for discovery of genes associated with detected sources of variation, gene annotation using publicly available databases and modules, and data visualization. 

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

* **Identify hidden sources of variability in transcriptomic data:**

![](https://github.com/nlawlor/iasva_shiny/blob/master/img/sv.plots.png)

* **Discover marker genes associated with a variable of interest:**

![](https://github.com/nlawlor/iasva_shiny/blob/master/img/marker.genes.png)

* **Determine molecular pathways associated with the identified marker genes:**

![](https://github.com/nlawlor/iasva_shiny/blob/master/img/pathway.analysis.png)

* **Interactive visualization of data:**

![](https://github.com/nlawlor/iasva_shiny/blob/master/img/tsne.gif)

***

### Availability

* **VIA-SVA** is currently hosted on **shinyapps.io** and available for use [here](https://nlawlor.shinyapps.io/VIASVA/).


### Other Resources

* **IA-SVA** is also available as an R-package that may be downloaded from **Bioconductor** [here](https://www.bioconductor.org/packages/devel/bioc/html/iasva.html), and the source code is available [here](https://github.com/UcarLab/iasva)

* The preprint manuscript describing the IA-SVA software can be found on **bioRxiv** [here](https://www.biorxiv.org/content/early/2018/04/24/151217)
