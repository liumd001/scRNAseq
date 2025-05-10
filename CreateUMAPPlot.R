########################################################
##    Create UMAP clustering Plots
#3    2023-02-21
########################################################

library(Seurat, quietly = TRUE)
library(tidyverse, quietly = TRUE)
require(DESeq2, quietly = TRUE)
require(GGally)
require(flextable)
require(ggplot2)
require(ggven)
require(VennDiagram)
require(patchwork)
require(Matrix)

out_dir = "U:\\team_share\\tfls\\2022_12_31_scRNASeq"

#tuong
tuong <- readRDS("W:/devtm/cbd/users/bdecato/Prostate Single Cell/PMID_34936871_fig1.rds")
tuong1 = tuong
tuong = tuong1

VlnPlot(tuong, features = c("nFeature_RNA","nCount_RNA"), ncol = 3)
plot1 = FeatureScatter(tuong, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")

#normalize
tuong = NormalizeData(tuong, normalization.method = "LogNormalize", scale.factor = 10000)

#identifying highly variable features
tuong = FindVariableFeatures(tuong, selection.method = "vst", nfeatures = 200)
top10 = head(VariableFeatures(tuong),10)
plot1 = VariableFeaturePlot(tuong)
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1  + plot2

#scaling data
all.genes = rownames(tuong)
tuong = ScaleData(tuong, features = all.genes)

#dimensional reduction
tuong = RunPCA(tuong, features = VariableFeatures(object = tuong))
print(tuong[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(tuong, dims = 1:2, reduction = "pca")
DimPlot(tuong, reduction = "pca")
DimHeatmap(tuong, dims = 1:15, cells = 500, balanced = TRUE)

#determine dimensionality
tuong = JackStraw(tuong, num.replicate = 100)
tuong = ScoreJackStraw(tuong, dims = 1:20)
JackStrawPlot(tuong, dims = 1:15)
ElbowPlot(tuong)

#clustering
tuong = FindNeighbors(tuong, dims = 1:10)
tuong = FindClusters(tuong, resolution = 0.5)
head(Idents(tuong),5)

#UMAP
tuong = RunUMAP(tuong, dims = 1:14)
DimPlot(tuong, reduction = "umap", group.by = 'Author.s.celltyp0e')
DimPlot(tuong, label = T)

g_umap_tuong  = tuong[['umap']]@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column(var = "cellid") %>%
  left_join(tuong@meta.data %>%
              select(Author.s.cell.type) %>%
              rownames_to_column(var = "cellid"), by = "cellid") %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = Author.s.cell.type)) +
  geom_point() +
  theme_bw() 

VlnPlot(tuong, features = c("KLK3","STEAP1"))
FeaturePlot(tuong, features = c("KLK3", "GAPDH", "FOLH1", "STEAP1", "STEAP2"))

#chen
chen <- readRDS("W:/devtm/cbd/users/bdecato/Prostate Single Cell/GSE141445.rds")

#normalize
chen = NormalizeData(chen, normalization.method = "LogNormalize", scale.factor = 10000)

#identifying highly variable features
chen = FindVariableFeatures(chen, selection.method = "vst", nfeatures = 200)
top10 = head(VariableFeatures(chen),10)
plot1 = VariableFeaturePlot(chen)
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1  + plot2

#scaling data
all.genes = rownames(chen)
chen = ScaleData(chen, features = all.genes)

#dimensional reduction
chen = RunPCA(chen, features = VariableFeatures(object = chen))
print(chen[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(chen, dims = 1:2, reduction = "pca")
DimPlot(chen, reduction = "pca")
DimHeatmap(chen, dims = 1:15, cells = 500, balanced = TRUE)

#determine dimensionality
chen = JackStraw(chen, num.replicate = 100)
chen = ScoreJackStraw(chen, dims = 1:20)
JackStrawPlot(chen, dims = 1:15)
ElbowPlot(chen)

#clustering
chen = FindNeighbors(chen, dims = 1:10)
chen = FindClusters(chen, resolution = 0.5)
head(Idents(chen),5)

#UMAP
chen = RunUMAP(chen, dims = 1:14)
DimPlot(chen, reduction = "umap", group.by = 'Author.s.cell.type')

g_umap_chen  = chen[['umap']]@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column(var = "cellid") %>%
  left_join(chen@meta.data %>%
              select(Author.s.cell.type) %>%
              rownames_to_column(var = "cellid"), by = "cellid") %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = Author.s.cell.type)) +
  geom_point() +
  theme_bw() 

#Song
dge = readRDS("C:\\Users\\mliu10\\Documents\\Amgen_documents\\prostate_scRNAseq\\GSE176031\\Seurat_dataset\\dge_E.rds")
DimPlot(dge, reduction = "umap")
g_umap_song  = dge[['umap']]@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column(var = "cellid") %>%
  left_join(dge@meta.data %>%
              select(ID) %>%
              rownames_to_column(var = "cellid"), by = "cellid") %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = ID)) +
  labs(color = "Cell Type") +
  geom_point() +
  theme_bw() 

