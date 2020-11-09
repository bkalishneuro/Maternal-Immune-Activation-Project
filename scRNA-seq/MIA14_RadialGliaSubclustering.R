25#### Seurat Main Cluster analysis on Orchestra
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)

#Inputs
species="m"
base_dir <- "/n/scratch2/bf78/MIASingleCell/MIA14REDO"
projectName = "MIA14REDO_BC_RG"
mapping_folder <- "/n/scratch2/bf78/MIASingleCell/MIA14REDO"

#Inputs
projectName = "MIA14REDO_BC_RG"
pc_num = 10;
res=0.4;
min_genes = 500;
load("cd_s_comb.Robj")

#Create a vector of Radial Glia Cells by selecting only clusters with radial glia marker gene Vim expressed
cd_comb_RG <- cd_comb[,RGCells]

seurat_RG <- CreateSeuratObject(counts = cd_comb_RG, min.cells = 3, project = projectName)

seurat_RG[["percent.mt"]] <- PercentageFeatureSet(seurat_RG, pattern = "^MT-")

head(seurat_RG@meta.data, 5)

seurat_RG <- subset(seurat_RG, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#Normalize
seurat_RG <- NormalizeData(object = seurat_RG, normalization.method = "LogNormalize", scale.factor = 10000)
print("Normalization is completed")

#ID variable genes
seurat_RG <- FindVariableFeatures(seurat_RG, selection.method = "vst", nfeatures = 2000)
print("FindVariableGenes is completed")

#Regression and Scaling with 8 cores!
all.genes <- rownames(seurat_RG)
seurat_RG <- ScaleData(seurat_RG, features = all.genes)
print("ScaleData and Regression is completed")

#Run PCA
seurat_RG <- RunPCA(seurat_RG, features = VariableFeatures(object = seurat_RG))

pdf("Elbowplot.pdf")
ElbowPlot(seurat_RG)
dev.off()

##Find clusters, VlnPlot per cluster, and see whether some of your clusters are formed of cells that systematically have lower or higher number of expressed genes
seurat_RG <- FindNeighbors(seurat_RG, dims = 1:10)
seurat_RG <- FindClusters(seurat_RG, resolution = res)


#Run tSNE and save coordinates
seurat_RG <- RunUMAP(seurat_RG, dims = 1:10,min.dist = 0.01)

embeds = Embeddings(seurat_RG[["umap"]])

write.csv(embeds, file = paste("seurat_RG_umapcoordinates",10,"_",res,".csv",sep=""))

seuratClusters <- Idents(seurat_RG)

write.csv(seuratClusters, file = "seurat_RG_seuratClusters.csv")

print("RunUmap is done")


ClusterMarkers <- FindAllMarkers(seurat_RG, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ClusterMarkers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(ClusterMarkers, file = "Final_ClusterMarkers_RG.csv")

save(seurat_RG, file = "seurat_RG_E14REDO.Robj")

print("seurat_RG has been saved")


