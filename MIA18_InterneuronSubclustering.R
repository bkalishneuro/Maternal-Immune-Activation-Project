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
res=0.3;
min_genes = 500;
load("seurat_mat_NewFinalE18.Robj")
load("cd_s_comb.Robj")
#Create a new folder
#Create CellVector, vector of all cell names listed as character strings, in this case all cells that express interneuron marker Gad2
SISeurat <- subset(seurat_mat_NewFinalE18,idents = c("1","3","4","7","9","11"))
SICells <- colnames(SISeurat)
NewComb <- cd_comb[,SICells]

### Seurat
Seurat_Striatal_Interneurons <- CreateSeuratObject(counts = NewComb, min.cells = 3, project = projectName)

Seurat_Striatal_Interneurons[["percent.mt"]] <- PercentageFeatureSet(Seurat_Striatal_Interneurons, pattern = "^MT-")

Seurat_Striatal_Interneurons[["percent.ribo"]] <- PercentageFeatureSet(Seurat_Striatal_Interneurons, pattern = "^RP[SL]")

head(Seurat_Striatal_Interneurons@meta.data, 5)

Seurat_Striatal_Interneurons <- subset(Seurat_Striatal_Interneurons, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#Normalize
Seurat_Striatal_Interneurons <- NormalizeData(object = Seurat_Striatal_Interneurons, normalization.method = "LogNormalize", scale.factor = 10000)
print("Normalization is completed")

#ID variable genes
Seurat_Striatal_Interneurons <- FindVariableFeatures(Seurat_Striatal_Interneurons, selection.method = "vst", nfeatures = 2000)
print("FindVariableGenes is completed")

#Regression and Scaling with 8 cores!
Seurat_Striatal_Interneurons <- ScaleData(Seurat_Striatal_Interneurons,vars.to.regress = c("percent.ribo","percent.mt"))
print("ScaleData and Regression is completed")

#Run PCA
Seurat_Striatal_Interneurons <- RunPCA(Seurat_Striatal_Interneurons, features = VariableFeatures(object = Seurat_Striatal_Interneurons))

pdf("Elbowplot.pdf")
ElbowPlot(Seurat_Striatal_Interneurons)
dev.off()

##Find clusters, VlnPlot per cluster, and see whether some of your clusters are formed of cells that systematically have lower or higher number of expressed genes
Seurat_Striatal_Interneurons <- FindNeighbors(Seurat_Striatal_Interneurons, dims = 1:20)
Seurat_Striatal_Interneurons <- FindClusters(Seurat_Striatal_Interneurons, resolution = res)

save(Seurat_Striatal_Interneurons, file = "Seurat_Striatal_Interneurons.Robj")

saveRDS(Seurat_Striatal_Interneurons, file = "Seurat_Striatal_Interneurons.rds")

#Run UMAP and save coordinates
Seurat_Striatal_Interneurons <- RunUMAP(Seurat_Striatal_Interneurons, dims = 1:20, min.dist = 0.01)

embeds = Embeddings(Seurat_Striatal_Interneurons[["umap"]])

write.csv(embeds, file = paste("umapcoordinatesE18StriatalInterneurons",30,"_",res,".csv",sep=""))

seuratClusters <- Idents(Seurat_Striatal_Interneurons)

write.csv(seuratClusters, file = "seuratClusterE18StriatalInterneurons.csv")

pdf("umapE18Striatalinterneurons.pdf")
DimPlot(Seurat_Striatal_Interneurons, reduction = "umap", label = TRUE, pt.size = 0.3) + NoLegend()
dev.off()
print("RunUmap is done")


ClusterMarkers <- FindAllMarkers(Seurat_Striatal_Interneurons, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ClusterMarkers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(ClusterMarkers, file = "ClusterMarkersE18Striatal_Interneurons.csv")






