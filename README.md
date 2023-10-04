# Cirrhotic lipid project: single-cell gene expression in liver
A GitHub containing the code for analysing single-cell expression of candidate genes (Cebola lab + Zoe Hall collaboration)

Corresponding author: Hannah Maude, hannah.maude12@imperial.ac.uk

**Plots to create:**

A) [Expression by cell type](#expression-by-cell-type): expression per cell population dot plot of set of 30 genes (Guilliams et al. dataset)

B) [Candidate gene UMAPs](#candidate-gene-UMAPs): some individual UMAPS of key genes Guilliams data set)

C) [Cirrhotic vs healthy expression](#cirrhotic-vs-healthy-expression): some cirrhotic v healthy violin plots from Ramachandran dataset of a few key genes and their relevant cell populations.

## Expression by cell type

First, we will explore the expression of our candidate genes using the scRNA-seq data from [Guilliams et al. (2022)](https://www.cell.com/cell/fulltext/S0092-8674(21)01481-1), available via the [Liver Cell Atlas](https://www.livercellatlas.org/datasets_human.php). 

The genes of interest:
| Gene  | Function |
| ----- | ---- |
| LPCAT1 | Lipid remodelling |
| LPCAT2 | Lipid remodelling |
FAR1	| Plasmalogen synthesis
SPTLC2	| Ceramide metabolism
CERS5	| Ceramide metabolism
CERS6	| Ceramide metabolism
DEGS1	| Ceramide metabolism
DEGS2	| Ceramide metabolism
CERK	| Ceramide metabolism
SPHK1	| Ceramide metabolism
SGPP2	| Ceramide metabolism
SGMS1	| Ceramide metabolism
SMPD1	| Ceramide metabolism
SMPD3	| Ceramide metabolism
ASAH1	| Ceramide metabolism
UGCG	| Glycosphingolipid metabolism
GBA	| Glycosphingolipid metabolism
B4GALT5	| Glycosphingolipid metabolism
B4GALT6	| Glycosphingolipid metabolism
ARSA	| Glycosphingolipid metabolism
A4GALT	| Glycosphingolipid metabolism
GBGT1	| Glycosphingolipid metabolism
B3GALT4	| Glycosphingolipid metabolism
B3GNT5	| Glycosphingolipid metabolism
GLA	| Glycosphingolipid metabolism
HEXA	| Glycosphingolipid metabolism
NEU1    | Glycosphingolipid metabolism

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

| ID | orig.ident | nCount_RNA | nFeature_RNA | cluster| annot | nCount_SCT | nFeature_SCT |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | 
| AAACCTGTCAGGATCT-1 |	SeuratProject | 4080 | 1343 | 34 | Mono+mono  derived cells | 3528 | 1343 | 
AAACGGGAGCACAGGT-1 | SeuratProject	| 7833	| 1857 |	34 |	Mono+mono derived cells	| 3648	| 698 |
AAACGGGAGTGTGGCA-1 | SeuratProject | 2665 | 1019 | 14 | Mono+mono derived cells	| 2724 | 1019 |
AAACGGGAGTTAGGTA-1 | SeuratProject | 1615 | 719	| 7 | Mono+mono derived cells	| 2402 | 719 |
AAACGGGCACCAGATT-1 | SeuratProject | 2537 | 848	| 7 | Mono+mono derived cells | 2652 | 848 |
AAACGGGTCATTGCCC-1 | SeuratProject | 5905 | 1406 | 43 | cDC1s | 4157	| 1403 | 

```r
# extract the meta data and add the UMAP coordinates
umapCoord = cell.annot[,1:2]
rownames(umapCoord) = cell.annot$cell
umapCoord = umapCoord[rownames(liver.filtered@meta.data),]
liver.filtered = AddMetaData(liver.filtered, metadata = umapCoord)
```

```r
# add it by CreateDimReducObject
liver.filtered[['umap']] <- CreateDimReducObject(embeddings = as.matrix(liver.filtered@meta.data[,c('UMAP_1','UMAP_2')]),
                                key = 'umap_', assay = 'RNA')
DimPlot(object = liver.filtered, label = TRUE,group.by = 'annot', reduction = "umap") + NoLegend() #+ ggtitle("Integrated controls") #+ ylim(-12,15) + xlim(-20,10)
```

<img src="https://github.com/CebolaLab/Cirrhotic_lipids/blob/main/Figures/Liver_atlas_all_clusters_github.png" width="50%" height="50%">


Read in the file containing the genes of interest:

```r
genes = read.delim('candidate_genes_dotplots.txt')
head(genes)
```
| Gene | Function |
| ---- | ---- | 
ELOVL1	| Lipid remodelling
LPCAT1	| Lipid remodelling
LPCAT2	| Lipid remodelling
FAR1	| Lipid remodelling
GNPAT	| Lipid remodelling
SPTLC2	| Ceramide metabolism

Dotplot of expression by cell type:

```r
DotPlot(object = liver.filtered, features = unique(genes$Gene), group.by='annot') + theme(axis.text.x = element_text(angle = 45,hjust=1))
ggsave('genes_of_interest_dotplot.pdf',width=10)
#ggsave('genes_of_interest_dotplot.png',width=10)
```

<img src="https://github.com/CebolaLab/Cirrhotic_lipids/blob/main/Figures/genes_of_interest_dotplot.png" width="100%" height="100%">


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

```r
# Plot UMAPs using the FeaturePlot command
FeaturePlot(object = liver.filtered, features = 'PSAP', label=FALSE, max.cutoff = 2, slot = 'data)
ggsave('UMAP_PSAP_maxcutoff_2_small.pdf',width=3.5,height=3)
```
<img src="https://github.com/CebolaLab/Cirrhotic_lipids/blob/main/Figures/UMAP_PSAP_maxcutoff_4.png" width="50%" height="50%">

## Cirrhotic vs healthy expression