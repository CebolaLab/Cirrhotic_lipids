# Cirrhotic lipid project: single-cell gene expression in liver
A GitHub containing the code for analysing single-cell expression of candidate genes (Cebola lab + Zoe Hall collaboration)

Corresponding author: Hannah Maude, hannah.maude12@imperial.ac.uk

**Plots to create:**

A) [Expression by cell type](#expression-by-cell-type): expression per cell population dot plot of set of 30 genes (Guilliams et al. dataset)

B) [Candidate gene UMAPs](#candidate-gene-UMAPs): some individual UMAPS of key genes Guilliams data set)

C) some cirrhotic v healthy violin plots from Ramachandran dataset of a few key genes and their relevant cell populations.

## Expression by cell type

First, we will explore the expression of our candidate genes using the scRNA-seq data from [Guilliams et al. (2022)](https://www.cell.com/cell/fulltext/S0092-8674(21)01481-1), available via the [Liver Cell Atlas](https://www.livercellatlas.org/datasets_human.php). 

The genes of interest:
| Gene  | Function |
| ----- | ---- |
| ELOVL1 | Lipid remodelling |
| LPCAT1 | Lipid remodelling |
| LPCAT2 | Lipid remodelling |
FAR1	| Lipid remodelling
GNPAT	| Lipid remodelling
SPTLC2	| Ceramide metabolism
CERS2	| Ceramide metabolism
CERS5	| Ceramide metabolism
CERS6	| Ceramide metabolism
DEGS1	| Ceramide metabolism
CERK	| Ceramide metabolism
SPHK1	| Ceramide metabolism
SGPP2	| Ceramide metabolism
SGMS1	| Ceramide metabolism
SMPD1	| Ceramide metabolism
SMPD3	| Ceramide metabolism
ASAH1	| Ceramide metabolism
S1PR4	| Signalling
UGCG	| Glycosphingolipid metabolism
PSAP	| Glycosphingolipid metabolism
B4GALT5	| Glycosphingolipid metabolism
B4GALT6	| Glycosphingolipid metabolism
ARSA	| Glycosphingolipid metabolism
A4GALT	| Glycosphingolipid metabolism
GBGT1	| Glycosphingolipid metabolism
B3GALT4	| Glycosphingolipid metabolism
B3GNT5	| Glycosphingolipid metabolism
GLA	| Glycosphingolipid metabolism
GBA	| Glycosphingolipid metabolism
HEXA	| Glycosphingolipid metabolism

```r
# Load libraries
library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')
library('limma')
library('ggplot2')
```

The single-cell RNA-seq data from the [Guilliams et al. 2022](https://www.cell.com/cell/fulltext/S0092-8674(21)01481-1) Liver Cell Atlas is available at this link: [https://www.livercellatlas.org/datasets_human.php](https://www.livercellatlas.org/datasets_human.php). The data can be downloaded from [this page](https://www.livercellatlas.org/download.php), including the gene-cell count matrix and cell annotation matrix for all liver cells, or analysis subsets (myeloid cells, lymphoid cells, CD45- cells). Download the *cell annotation matrix* for all liver cells:

```bash
# Run this on the command line (bash)
wget https://www.livercellatlas.org/data_files/toDownload/annot_humanAll.csv
```

```r
# The data can then be read into R
# Use the read.table command to read in the annotation file from the current directory and save it as a dataframe called cell.annot
cell.annot = read.table('annot_humanAll.csv',sep=',',header=TRUE)
```

> NOTE: the gene-cell count matrices contain all genes and all cells before any quality control filtering was done. The cell annotation matrices contain the cells that were retained after QC filtering.

Download the cell count matrices and extract the zipped file (rawData_human.zip). 

```bash 
# This is run on the command line (bash)
wget https://www.livercellatlas.org/data_files/toDownload/rawData_human.zip
unzip rawData_human.zip
```

The data can be read into R:
```r
# Read in the count tables containing the raw data
liverall = Read10X('rawData_human/countTable_human/',
  gene.column = 1,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE)
  
# Create a Seurat object using the data
liverall_object = CreateSeuratObject(counts = liverall)
```

Filter the `liverall` object using the `cell.annot` dataframe, which contains cell retained after QC filtering.

```r
liver.filtered = subset(liverall_object, cells = cell.annot$cell)
# Add information to the metadata in the liver.filtered object.
liver.filtered[['cluster']] = cell.annot$cluster
liver.filtered[['annot']] = cell.annot$annot
```

Normalise the raw data using the Seurat `SCTransform` command.

```r
#Look up the SCTransform on the Seurat website. 
liver.filtered = SCTransform(liver.filtered, conserve.memory = TRUE,return.only.var.genes = FALSE) 
```

```r
head(liver.filtered@meta.data)
```

| | orig.ident | nCount_RNA | nFeature_RNA | cluster| annot | nCount_SCT | nFeature_SCT |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | 
| AAACCTGTCAGGATCT-1 |	SeuratProject | 4080 | 1343 | 34 | Mono+mono  derived cells | 3528 | 1343 | 
AAACGGGAGCACAGGT-1 | SeuratProject	| 7833	| 1857 |	34 |	Mono+mono derived cells	| 3648	| 698 |
AAACGGGAGTGTGGCA-1 | SeuratProject | 2665 | 1019 | 14 | Mono+mono derived cells	| 2724 | 1019 |
AAACGGGAGTTAGGTA-1 | SeuratProject | 1615 | 719	| 7 | Mono+mono derived cells	| 2402 | 719 |
AAACGGGCACCAGATT-1 | SeuratProject | 2537 | 848	| 7 | Mono+mono derived cells | 2652 | 848 |
AAACGGGTCATTGCCC-1 | SeuratProject | 5905 | 1406 | 43 | cDC1s | 4157	| 1403 | 



## Candidate gene UMAPs
UMAP:
| Gene	| B |
| ---- | ---- | 
LPCAT1 |Endothelial and NKT
LPCAT2	| Macrophages
FAR1	| Endothelial, stromal, other
SPTLC2	| Variety - highest in macropahge, endothelial
SGMS1	| Hep, stromal and others
ASAH1	| Macrophage and others
S1PR4	| Monocytes, NK, T
UGCG	| Dendirtic and others
PSAP	| High esp mono and macro