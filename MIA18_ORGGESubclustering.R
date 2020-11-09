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

#Create CellVector, vector of all cell names listed as character strings, in this case all clusters that express markers for oligodendrocytes, radial glia, or ganglionic eminence
ORGGESeurat <- subset(seurat_mat_NewFinalE18,idents = c("12","13","14","15","19"))
ORGGECells <- colnames(ORGGESeurat)
NewComb <- cd_comb[,ORGGECells]

### Seurat
Seurat_Oligo_RG_GE <- CreateSeuratObject(counts = NewComb, min.cells = 3, project = projectName)

Seurat_Oligo_RG_GE[["percent.mt"]] <- PercentageFeatureSet(Seurat_Oligo_RG_GE, pattern = "^MT-")

Seurat_Oligo_RG_GE[["percent.ribo"]] <- PercentageFeatureSet(Seurat_Oligo_RG_GE, pattern = "^RP[SL]")

head(Seurat_Oligo_RG_GE@meta.data, 5)

Seurat_Oligo_RG_GE <- subset(Seurat_Oligo_RG_GE, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#Normalize
Seurat_Oligo_RG_GE <- NormalizeData(object = Seurat_Oligo_RG_GE, normalization.method = "LogNormalize", scale.factor = 10000)
print("Normalization is completed")

#ID variable genes
Seurat_Oligo_RG_GE <- FindVariableFeatures(Seurat_Oligo_RG_GE, selection.method = "vst", nfeatures = 2000)
print("FindVariableGenes is completed")

#Regression and Scaling with 8 cores!
Seurat_Oligo_RG_GE <- ScaleData(Seurat_Oligo_RG_GE,vars.to.regress = c("percent.ribo","percent.mt"))
print("ScaleData and Regression is completed")

#Run PCA
Seurat_Oligo_RG_GE <- RunPCA(Seurat_Oligo_RG_GE, features = VariableFeatures(object = Seurat_Oligo_RG_GE))

pdf("Elbowplot.pdf")
ElbowPlot(Seurat_Oligo_RG_GE)
dev.off()

##Find clusters, VlnPlot per cluster, and see whether some of your clusters are formed of cells that systematically have lower or higher number of expressed genes
Seurat_Oligo_RG_GE <- FindNeighbors(Seurat_Oligo_RG_GE, dims = 1:20)
Seurat_Oligo_RG_GE <- FindClusters(Seurat_Oligo_RG_GE, resolution = res)

save(Seurat_Oligo_RG_GE, file = "Seurat_Oligo_RG_GE.Robj")

saveRDS(Seurat_Oligo_RG_GE, file = "Seurat_Oligo_RG_GE.rds")

#Run UMAP and save coordinates
Seurat_Oligo_RG_GE <- RunUMAP(Seurat_Oligo_RG_GE, dims = 1:20, min.dist = 0.01)

embeds = Embeddings(Seurat_Oligo_RG_GE[["umap"]])

write.csv(embeds, file = paste("umapcoordinatesE18Oligo_RG_GE",30,"_",res,".csv",sep=""))

seuratClusters <- Idents(Seurat_Oligo_RG_GE)

write.csv(seuratClusters, file = "seuratClusterE18Oligo_RG_GE.csv")

pdf("umapE18Oligo_RG_GE.pdf")
DimPlot(Seurat_Oligo_RG_GE, reduction = "umap", label = TRUE, pt.size = 0.3) + NoLegend()
dev.off()
print("RunUmap is done")


ClusterMarkers <- FindAllMarkers(Seurat_Oligo_RG_GE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ClusterMarkers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(ClusterMarkers, file = "ClusterMarkersE18Oligo_RG_GE.csv")






