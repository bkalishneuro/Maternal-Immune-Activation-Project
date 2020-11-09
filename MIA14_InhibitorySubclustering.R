### Seurat Main Cluster analysis on O2, E14 Int Subclustering

library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
projectName = "MIA14REDO_BC_Interneuron"
pc_num = 30;
res=0.6;
min_genes = 500;


#Create a vector of Int Cells by selecting only clusters with Int marker gene Gad2
cd_comb_interneurons <- cd_comb[,Interneurons]

seurat_Int <- CreateSeuratObject(counts = cd_comb_Int, min.cells = 3, project = projectName)

seurat_Int[["percent.mt"]] <- PercentageFeatureSet(seurat_Int, pattern = "^MT-")

head(seurat_Int@meta.data, 5)

seurat_Int <- subset(seurat_Int, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#Normalize
seurat_Int <- NormalizeData(object = seurat_Int, normalization.method = "LogNormalize", scale.factor = 10000)
print("Normalization is completed")

#ID variable genes
seurat_Int <- FindVariableFeatures(seurat_Int, selection.method = "vst", nfeatures = 2000)
print("FindVariableGenes is completed")

#Regression and Scaling with 8 cores!
all.genes <- rownames(seurat_Int)
seurat_Int <- ScaleData(seurat_Int, features = all.genes)
print("ScaleData and Regression is completed")

#Run PCA
seurat_Int <- RunPCA(seurat_Int, features = VariableFeatures(object = seurat_Int))

pdf("Elbowplot.pdf")
ElbowPlot(seurat_Int)
dev.off()

##Find clusters, VlnPlot per cluster, and see whether some of your clusters are formed of cells that systematically have lower or higher number of expressed genes
seurat_Int <- FindNeighbors(seurat_Int, dims = 1:30)
seurat_Int <- FindClusters(seurat_Int, resolution = res)

#Run UMAP and save coordinates
seurat_Int <- RunUMAP(seurat_Int, dims = 1:30,min.dist = 0.1)

embeds = Embeddings(seurat_Int[["umap"]])

write.csv(embeds, file = paste("seurat_Int_umapcoordinates",30,"_",res,".csv",sep=""))

seuratClusters <- Idents(seurat_Int)

write.csv(seuratClusters, file = "seurat_Int_seuratClusters.csv")

pdf("FullyFilteredE14umap.pdf")
DimPlot(seurat_Int, reduction = "umap",label = TRUE)
dev.off()
print("RunUmap is done")


ClusterMarkers <- FindAllMarkers(seurat_Int, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ClusterMarkers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(ClusterMarkers, file = "ClusterMarkersE14_Int.csv")

save(seurat_Int, file = "seurat_Int_E14REDO.Robj")

print("seurat_Int has been saved")


