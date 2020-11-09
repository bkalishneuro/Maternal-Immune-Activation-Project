### Seurat Main Cluster analysis on O2, E14 Cortical Subclustering

library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
projectName = "MIA14REDO_BC_Cortical"
pc_num = 30;
res=0.6;
min_genes = 500;

#Create a vector of Excitatory Cortical Cells by selecting only clusters with excitatory cortical cell marker genes Neurod2 and Neurod6
cd_comb_cortical <- cd_comb[,CorticalCells]

seurat_Cortical <- CreateSeuratObject(counts = cd_comb_cortical, min.cells = 3, project = projectName)

seurat_Cortical[["percent.mt"]] <- PercentageFeatureSet(seurat_Cortical, pattern = "^MT-")

head(seurat_Cortical@meta.data, 5)

seurat_Cortical <- subset(seurat_Cortical, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#Normalize
seurat_Cortical <- NormalizeData(object = seurat_Cortical, normalization.method = "LogNormalize", scale.factor = 10000)
print("Normalization is completed")

#ID variable genes
seurat_Cortical <- FindVariableFeatures(seurat_Cortical, selection.method = "vst", nfeatures = 2000)
print("FindVariableGenes is completed")

#Regression and Scaling with 8 cores!
all.genes <- rownames(seurat_Cortical)
seurat_Cortical <- ScaleData(seurat_Cortical, features = all.genes)
print("ScaleData and Regression is completed")

#Run PCA
seurat_Cortical <- RunPCA(seurat_Cortical, features = VariableFeatures(object = seurat_Cortical))

pdf("Elbowplot.pdf")
ElbowPlot(seurat_Cortical)
dev.off()

##Find clusters, VlnPlot per cluster, and see whether some of your clusters are formed of cells that systematically have lower or higher number of expressed genes
seurat_Cortical <- FindNeighbors(seurat_Cortical, dims = 1:30)
seurat_Cortical <- FindClusters(seurat_Cortical, resolution = res)


#Run UMAP and save coordinates
seurat_Cortical <- RunUMAP(seurat_Cortical, dims = 1:30,min.dist = 0.1)

embeds = Embeddings(seurat_Cortical[["umap"]])

write.csv(embeds, file = paste("seurat_Cortical_umapcoordinates",30,"_",res,".csv",sep=""))

seuratClusters <- Idents(seurat_Cortical)

write.csv(seuratClusters, file = "seurat_Cortical_seuratClusters.csv")

pdf("FullyFilteredE14umap.pdf")
DimPlot(seurat_Cortical, reduction = "umap",label = TRUE)
dev.off()
print("RunUmap is done")


ClusterMarkers <- FindAllMarkers(seurat_Cortical, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ClusterMarkers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(ClusterMarkers, file = "ClusterMarkersE14_Cortical.csv")

save(seurat_Cortical, file = "seurat_Cortical_E14REDO.Robj")

print("seurat_Cortical has been saved")


