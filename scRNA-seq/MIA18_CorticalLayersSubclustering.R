#### Seurat Main Cluster analysis on Orchestra
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)

#Inputs
species="m"
base_dir <- "/n/scratch2/bf78/MIASingleCell/MIA18"
projectName = "MIA18"

pc_num = 20;
res=0.4;
min_genes = 500;

load("seurat_mat_NewFinalE18.Robj")
load("cd_s_comb.Robj")

#Create CellVector, vector of all cell names listed as character strings, in this case, all cells expressing excitatory coritcal layer neuron marker genes Neurod2 and Neurod 6
CLSeurat <- subset(seurat_mat_NewFinalE18,idents = c("0","2","5","6","8","10","16","20"))
CLCells <- colnames(CLSeurat)
NewComb <- cd_comb[,CLCells]

### Seurat
Seurat_CorticalLayers <- CreateSeuratObject(counts = NewComb, min.cells = 3, project = projectName)

Seurat_CorticalLayers[["percent.mt"]] <- PercentageFeatureSet(Seurat_CorticalLayers, pattern = "^MT-")

Seurat_CorticalLayers[["percent.ribo"]] <- PercentageFeatureSet(Seurat_CorticalLayers, pattern = "^RP[SL]")

head(Seurat_CorticalLayers@meta.data, 5)

Seurat_CorticalLayers <- subset(Seurat_CorticalLayers, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#Normalize
Seurat_CorticalLayers <- NormalizeData(object = Seurat_CorticalLayers, normalization.method = "LogNormalize", scale.factor = 10000)
print("Normalization is completed")

#ID variable genes
Seurat_CorticalLayers <- FindVariableFeatures(Seurat_CorticalLayers, selection.method = "vst", nfeatures = 2000)
print("FindVariableGenes is completed")

#Regression and Scaling with 8 cores!
Seurat_CorticalLayers <- ScaleData(Seurat_CorticalLayers,vars.to.regress = c("percent.ribo","percent.mt"))
print("ScaleData and Regression is completed")

#Run PCA
Seurat_CorticalLayers <- RunPCA(Seurat_CorticalLayers, features = VariableFeatures(object = Seurat_CorticalLayers))

pdf("Elbowplot.pdf")
ElbowPlot(Seurat_CorticalLayers)
dev.off()

##Find clusters, VlnPlot per cluster, and see whether some of your clusters are formed of cells that systematically have lower or higher number of expressed genes
Seurat_CorticalLayers <- FindNeighbors(Seurat_CorticalLayers, dims = 1:20)
Seurat_CorticalLayers <- FindClusters(Seurat_CorticalLayers, resolution = res)

save(Seurat_CorticalLayers, file = "Seurat_CorticalLayers.Robj")

saveRDS(Seurat_CorticalLayers, file = "Seurat_CorticalLayers.rds")

#Run UMAP and save coordinates
Seurat_CorticalLayers <- RunUMAP(Seurat_CorticalLayers, dims = 1:20, min.dist = 0.01)

embeds = Embeddings(Seurat_CorticalLayers[["umap"]])

write.csv(embeds, file = paste("umapcoordinatesE18CorticalLayers",30,"_",res,".csv",sep=""))

seuratClusters <- Idents(Seurat_CorticalLayers)

write.csv(seuratClusters, file = "seuratClusterE18CorticalLayers.csv")

pdf("umapE18CorticalLayers.pdf")
DimPlot(Seurat_CorticalLayers, reduction = "umap", label = TRUE, pt.size = 0.3) + NoLegend()
dev.off()
print("RunUmap is done")


ClusterMarkers <- FindAllMarkers(Seurat_CorticalLayers, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ClusterMarkers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(ClusterMarkers, file = "ClusterMarkersE18CorticalLayers.csv")






