# Gruevska et al: single-cell gene expression in liver
A GitHub containing the code for analysing single-cell expression of candidate genes (Gruevska et al).

Corresponding author: Hannah Maude, hannah.maude12@imperial.ac.uk

**Plots to create:**

A) [Expression by cell type](#expression-by-cell-type): expression per cell population dot plot of set of 30 genes (Liver Cell Atlas: Guilliams et al. dataset)

B) [Candidate gene UMAPs](#candidate-gene-UMAPs): some individual UMAPS of key genes

C) [Cirrhotic vs healthy expression](#cirrhotic-vs-healthy-expression): cirrhotic vs healthy violin plots from the Ramachandran et al. dataset of a few key genes and their relevant cell populations.

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

The single-cell RNA-seq data from the [Guilliams et al. 2022](https://www.cell.com/cell/fulltext/S0092-8674(21)01481-1) Liver Cell Atlas is available at this link: [https://www.livercellatlas.org/datasets_human.php](https://www.livercellatlas.org/datasets_human.php). The data can be downloaded from [this page](https://www.livercellatlas.org/download.php), including the gene-cell count matrix and cell annotation matrix for all liver cells, or analysis subsets (myeloid cells, lymphoid cells, CD45- cells). The data, including the *cell annotation matrix*, can be downloaded for all liver cells via the command line:

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

The UMAP coordinates corresponding to the Liver Cell Atlas webpage can be added for consistency. Here, the UMAP will be recalculated to add the full model for reference-query dataset integration.

```r
# extract the meta data and add the UMAP coordinates
# umapCoord = cell.annot[,1:2]
# rownames(umapCoord) = cell.annot$cell
# umapCoord = umapCoord[rownames(liver.filtered@meta.data),]
# liver.filtered = AddMetaData(liver.filtered, metadata = umapCoord)
```

```r
liver.filtered <- RunPCA(liver.filtered, verbose = FALSE, return.model = TRUE)
liver.filtered <- RunUMAP(liver.filtered, reduction = "pca", dims = 1:30, return.model = TRUE)
```

```r
# add it by CreateDimReducObject
#liver.filtered[['umap']] <- CreateDimReducObject(embeddings = as.matrix(liver.filtered@meta.data[,c('UMAP_1','UMAP_2')]),key = 'umap_', assay = 'RNA')
DimPlot(object = liver.filtered, label = TRUE, group.by = 'annot', reduction = "umap") + NoLegend() #+ ggtitle("Integrated controls") #+ ylim(-12,15) + xlim(-20,10)
```

<img src="https://github.com/CebolaLab/Cirrhotic_lipids/blob/main/Figures/Liver_atlas_all_clusters_github.png" width="50%" height="50%">


Read in the file containing the genes of interest:

```r
genes = read.delim('main_gene_list_short.txt')
#genes = read.delim('main_gene_list_long.txt')
head(genes)
```
| Gene | Function |
| ---- | ---- | 
LPCAT1	| Lipid remodelling
LPCAT2	| Lipid remodelling
FAR1	| Lipid remodelling
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

Download Healthy and Cirrhotic sample files from GEO [GSE136103](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136103). Assuming these files are in the working directory. 

First, we will list the files in the directory with `_matrix.mtx.gz` in the name using `list.files()`, and create two vectors with the names of each sample by removing the `_matrix.mtx.gz` extension. 

```r
healthy.samples = gsub('_matrix.mtx.gz','',grep('cd45+',grep('healthy',grep('matrix',list.files(),value=TRUE),value=TRUE),value=TRUE))
cirrhotic.samples = gsub('_matrix.mtx.gz','',grep('cd45+',grep('cirrhotic',grep('matrix',list.files(),value=TRUE),value=TRUE),value=TRUE))
```

Next, we create a Seurat object for each sample using the `ReadMtx` file. This utilises a `for` loop, which runs the code for each sample in the vector of sample names created above. This command takes the `_matrix.mtx.gz`, `_barcodes.tsv.gz` and `_genes.tsv.gz` files as input. The code also replaces the + and - with "plus" or "neg".

```r
#Read data matrix files into R using the barcodes and gene dataframes provided. 
for(x in healthy.samples){
    data=ReadMtx(paste0(x,'_matrix.mtx.gz'),paste0(x,'_barcodes.tsv.gz'),paste0(x,'_genes.tsv.gz'))
    seurat_object = CreateSeuratObject(counts = data)
    if(grepl("+",x)){
        name = gsub('\\+',"plus",x)}
    if(grepl("-",x)){
        name = gsub('\\-',"neg",x)}
    assign(name, seurat_object)
}

for(x in cirrhotic.samples){
    data=ReadMtx(paste0(x,'_matrix.mtx.gz'),paste0(x,'_barcodes.tsv.gz'),paste0(x,'_genes.tsv.gz'))
    seurat_objec = CreateSeuratObject(counts = data)
    if(grepl("+",x)){
        name = gsub('\\+',"plus",x)}
    if(grepl("-",x)){
        name = gsub('\\-',"neg",x)}
    assign(name, seurat_object)
}

#Obtain a vector with the names of all the R objects you have created, which were in the list of sample names with '+' and '-' replaced with 'plus' and 'negative', respectively. 
all.samples = gsub('\\-','neg',gsub('\\+','plus',c(healthy.samples,cirrhotic.samples)))

#This applies the 'get' function to the vector of names to return the Seurat object of the same name, and saves them as a list.
all_list = lapply(all.samples,get)

#Assign the names of each item in the list as the sample names
names(all_list) = all.samples

#For each object (sample), add a metadata column with the name of the sample
all_list = mapply(function(all_list, i) AddMetaData(all_list,  metadata = as.character(i), col.name="orig.ident"), all_list, names(all_list),SIMPLIFY = FALSE)

#For each object (sample), add a metadata column with the phenotype (healthy or cirrhotic)
all_list = mapply(function(all_list, i) AddMetaData(all_list, metadata = substr(strsplit(as.character(i),'_')[[1]][2],1,nchar(strsplit(as.character(i),'_')[[1]][2])-1), col.name="phenotype") , all_list, names(all_list),SIMPLIFY = FALSE)   

library(glmGamPoi)
all_list <- lapply(X = all_list, FUN = SCTransform, vst.flavor = "v2")
# For each object (sample), run the PCA 
all_list = mapply(function(all_list) RunPCA(all_list,verbose = FALSE), all_list)
# For each object (sample), run the UMAP dimensionality reduction
all_list = mapply(function(all_list) RunUMAP(all_list, dims = 1:20, verbose = FALSE), all_list)
# For each object (sample), run the UMAP dimensionality reduction
all_list = mapply(function(all_list) FindNeighbors(all_list, dims = 1:20, verbose = FALSE), all_list)
# For each object (sample), run the UMAP dimensionality reduction
all_list = mapply(function(all_list) FindClusters(all_list, verbose = FALSE), all_list)
```

Detect and remove doublets.

```r
library(scDblFinder)

# Run scDblFinder for all the samples.
sce.list = list()
for(name in names(all_list)){
    sce.list[[name]] <- scDblFinder(as.SingleCellExperiment(all_list[[name]]), clusters="seurat_clusters")    
}
# Add the scores back to the Seurat object
for(name in names(all_list)){
    all_list[[name]]$scDblFinder.score <- sce.list[[name]]$scDblFinder.score
}
# Add the scores back to the Seurat object
for(name in names(all_list)){
    all_list[[name]]$scDblFinder.class <- sce.list[[name]]$scDblFinder.class
}
# For each object, subset to keep only singlets
for(name in names(all_list)){
    all_list[[name]] <- subset(all_list[[name]], subset = scDblFinder.class == 'singlet')
}
```

Repeat the processing now that doublets have been removed.

```r
all_list <- lapply(X = all_list, FUN = SCTransform, vst.flavor = "v2")
# For each object (sample), run the PCA 
all_list = mapply(function(all_list) RunPCA(all_list,verbose = FALSE), all_list)
# For each object (sample), run the UMAP dimensionality reduction
all_list = mapply(function(all_list) RunUMAP(all_list, dims = 1:20, verbose = FALSE), all_list)
# For each object (sample), run the UMAP dimensionality reduction
all_list = mapply(function(all_list) FindNeighbors(all_list, dims = 1:20, verbose = FALSE), all_list)
# For each object (sample), run the UMAP dimensionality reduction
all_list = mapply(function(all_list) FindClusters(all_list, verbose = FALSE), all_list)

```

Next, we filter to remove clusters of low quality cells. For each sample, cluster marker genes will be calculated and the cluster will be removed if there is a mitochondrial gene in the top five markers. 

```r
#For each object (sample), add a metadata column with the phenotype (healthy or cirrhotic)
all_list = mapply(function(all_list) AddMetaData(all_list, metadata = PercentageFeatureSet(all_list, pattern = "^MT-"),col.name="percent.mt"), all_list, SIMPLIFY = FALSE) 

all_list.markers = mapply(function(all_list) FindAllMarkers(all_list, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25), all_list, SIMPLIFY = FALSE)  

library(dplyr)

# Subset to the top five marker genes
for (i in 1:length(all_list.markers)) {
  # Subset the current dataframe based on a condition (e.g., 'gene' contains 'MT-')
  all_list.markers[[i]] <- all_list.markers[[i]] %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)
}

# For each object, extract the cluster IDs with 'MT-' in the top five.
filtered_list = list()
# Loop through each data frame in the list
for (i in 1:length(all_list.markers)) {
  # Get the name of the current dataframe
  df_name <- names(all_list.markers)[i]
   # Filter rows based on 'gene' column
  filtered_rows <- all_list.markers[[i]][grepl('MT-', all_list.markers[[i]]$gene), 'cluster']
  # Add the filtered rows to the result list with the dataframe name
  filtered_list[[df_name]] <- unique(as.character(unlist(filtered_rows)))
}

# subset for objects with > 0 clusters reported
filtered_list = subset(filtered_list, lapply(filtered_list,length)>0)

# For each sample with a MT cluster, remove the cluster
for(sample in names(filtered_list)){
    all_list[[sample]] <- subset(all_list[[sample]], ident = filtered_list[[sample]], invert = TRUE)
}
```

Repeat the processing and clustering with the cleaned data.

```r
all_list <- lapply(X = all_list, FUN = SCTransform, vst.flavor = "v2")
# For each object (sample), run the PCA 
all_list = mapply(function(all_list) RunPCA(all_list,verbose = FALSE), all_list)
# For each object (sample), run the UMAP dimensionality reduction
all_list = mapply(function(all_list) RunUMAP(all_list, dims = 1:20, verbose = FALSE), all_list)
# For each object (sample), run the UMAP dimensionality reduction
all_list = mapply(function(all_list) FindNeighbors(all_list, dims = 1:20, verbose = FALSE), all_list)
# For each object (sample), run the UMAP dimensionality reduction
all_list = mapply(function(all_list) FindClusters(all_list, verbose = FALSE), all_list)
```

Integrate the samples.

```r
# Try with a reduced number of features and PCA dimensions
# nfeatures reduced to 2500 to avoid errors
features <- SelectIntegrationFeatures(object.list = all_list, nfeatures = 2500)
all_list <- PrepSCTIntegration(object.list = all_list, anchor.features = features)
# dims reduced to 25 as using 30 returned an error
liver.anchors <- FindIntegrationAnchors(all_list, normalization.method = "SCT", anchor.features = features, dims = 1:25)

# this command creates an 'integrated' data assay
all.integrated <- IntegrateData(anchorset = liver.anchors,  normalization.method = "SCT", preserve.order = TRUE)

all.integrated <- RunPCA(all.integrated, verbose = FALSE)
all.integrated <- RunUMAP(all.integrated, reduction = "pca", dims = 1:30)

# specify that we will perform downstream analysis on the corrected data note that the original unmodified data still resides in the 'RNA' assay
DefaultAssay(all.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
all.integrated <- FindNeighbors(all.integrated, dims = 1:10)
all.integrated <- FindClusters(object = all.integrated) #, graph.name = NULL, algorithm = 4)

DimPlot(all.integrated, group.by = "orig.ident")
DimPlot(all.integrated, group.by = "phenotype")
DimPlot(all.integrated, group.by = "seurat_clusters", label=TRUE)
```

<img src="https://github.com/CebolaLab/Cirrhotic_lipids/blob/main/Figures/Ramachandran-umap1.png">

## Integrate the Ramachandran et al. dataset with the Liver Cell Atlas

```r
anchors <- FindTransferAnchors(reference = liver.filtered, query = all.integrated,
    dims = 1:30, reference.reduction = "pca")
    
all.integrated <- MapQuery(anchorset = anchors, reference = liver.filtered, query = all.integrated,
    refdata = list(celltype = "annot"), reference.reduction = "pca", reduction.model = "umap")
    
DimPlot(all.integrated, reduction = "ref.umap",  group.by = "predicted.celltype", split.by = 'phenotype', label = TRUE,label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
```

<img src="https://github.com/CebolaLab/Cirrhotic_lipids/blob/main/Figures/NEWest_Ramachandran_UMAP.png">

```r
all.integrated$celltype.phenotype <- paste(all.integrated$predicted.celltype, all.integrated$phenotype, sep = "_")

Idents(all.integrated) = 'predicted.celltype'

# For visualisation and differential gene expression, set the default assay back to SCT
DefaultAssay(all.integrated) <- "SCT"

# subset the query gene list for genes in the dataset
longgenes = subset(genes.long[,1], genes.long[,1] %in% rownames(all.integrated@assays$SCT@data))

dge = list()
for(x in unique(Idents(all.integrated))){
    tmp <- subset(all.integrated, idents = x, features = longgenes)
    Idents(tmp) <- 'phenotype'
    DefaultAssay(tmp) <- 'RNA'
    dge[[x]] <- FindMarkers(tmp, ident.1 = "cirrhotic", ident.2 = "healthy", verbose = FALSE)
}

# Define a function to add rownames as a 'Gene' column
add_gene_column <- function(df) {
  df$Gene <- rownames(df)
  return(df)
}

# Apply the function to each dataframe in df_list using mapply
dge <- mapply(add_gene_column, dge, SIMPLIFY = FALSE)

# Combine the dataframes and add a column with names
library(dplyr)
dge.combined = bind_rows(dge, .id = "celltype")

dge.combined[order(dge.combined$p_val_adj),]
```

<img src="https://github.com/CebolaLab/Cirrhotic_lipids/blob/main/Figures/degs_screenshot.png">


```r
celltype = c('Mono+mono derived cells','Neutrophils','Fibroblasts','Endothelial cells','Basophils','Mig.cDCs')
gene = 'ASAH1' #c('HEXA','HEXB','GLA')
tmp <- subset(all.integrated, idents = celltype) #,features="HEXA")
VlnPlot(tmp, features = gene, split.by = "phenotype", slot='data', split.plot = TRUE)
#ggsave(paste0(celltype,'_NEWviolin_',gene,'.pdf'))
ggsave('ASAH1_all_DEGs_NEWviolin.pdf')
```

<img src="https://github.com/CebolaLab/Cirrhotic_lipids/blob/main/Figures/ASAH1_all_DEGs_NEWviolin.png">