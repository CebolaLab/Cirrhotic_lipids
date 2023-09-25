# Cirrhotic_lipids
A GitHub containing the code for analysing single-cell expression of candidate genes (Cebola lab + Zoe Hall collaboration)

Corresponding author: Hannah Maude, hannah.maude12@imperial.ac.uk

**Plots to create:**

A) [Expression by cell type](#expression-by-cell-type): expression per cell population dot plot of set of 30 genes (Guilliams et al. dataset)

B) [Candidate gene UMAPs](#candidate-gene-UMAPs): some individual UMAPS of key genes Guilliams data set)

C) some cirrhotic v healthy violin plots from Ramachandran dataset of a few key genes and their relevant cell populations.

## Expression by cell type

Dot plot:
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

The single-cell RNA-seq data from the [Guilliams et al. 2022](https://www.cell.com/cell/fulltext/S0092-8674(21)01481-1) Liver Cell Atlas is available at this link: [https://www.livercellatlas.org/datasets_human.php](https://www.livercellatlas.org/datasets_human.php)

```r
#read.table is the command for reading a file in the current directory (folder) 
#cell.annot= means "save the file you have just read as a dataframe called cell.annot"
cell.annot=read.table('annot_humanAll.csv',sep=',',header=TRUE)
```

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