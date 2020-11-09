#### Seurat Main Cluster analysis on O2
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)

#Inputs
species="m"
base_dir <- "/n/groups/neuroduo/brian/MIAraw/"
projectName = "MIA18"
mapping_folder <- "/n/groups/neuroduo/brian/MIAraw/"

pc_num = 30;
res=1;
min_genes = 500;


########Load files
folder_names <- c("PBSME18")
if (species == "h") cd <- data.frame(test=c(1:25463)) else cd <- data.frame(test=c(1:25289));

dim(cd)

for (i in folder_names)
{
    counts <- read.delim(paste(mapping_folder,i,"/",i,".counts.tsv",sep=""),skip=0,header=T,sep="\t",stringsAsFactors=F,row.names=1,check.names = F);
    rownames(counts) <- paste(i,rownames(counts),sep = "_")
    counts <- counts[rowSums(counts)>min_genes,];
    counts_t <- as.data.frame(t(counts))
    cd <- cbind(cd,counts_t);
    print(i)
}

cd$test <- NULL
cd_s_1_1 <- as(as.matrix(cd), "dgCMatrix")
cd <- NULL
print("PBSME18")

#############Load files
folder_names <- c("PBS2ME18")
if (species == "h") cd <- data.frame(test=c(1:25463)) else cd <- data.frame(test=c(1:25289));

dim(cd)

for (i in folder_names)
{
    counts <- read.delim(paste(mapping_folder,i,"/",i,".counts.tsv",sep=""),skip=0,header=T,sep="\t",stringsAsFactors=F,row.names=1,check.names = F);
    rownames(counts) <- paste(i,rownames(counts),sep = "_")
    counts <- counts[rowSums(counts)>min_genes,];
    counts_t <- as.data.frame(t(counts))
    cd <- cbind(cd,counts_t);
    print(i)
}

cd$test <- NULL
cd_s_1_2 <- as(as.matrix(cd), "dgCMatrix")
cd <- NULL
print("PBS2ME18")

########Load files
folder_names <- c("PBSFE18")
if (species == "h") cd <- data.frame(test=c(1:25463)) else cd <- data.frame(test=c(1:25289));

dim(cd)

for (i in folder_names)
{
    counts <- read.delim(paste(mapping_folder,i,"/",i,".counts.tsv",sep=""),skip=0,header=T,sep="\t",stringsAsFactors=F,row.names=1,check.names = F);
    rownames(counts) <- paste(i,rownames(counts),sep = "_")
    counts <- counts[rowSums(counts)>min_genes,];
    counts_t <- as.data.frame(t(counts))
    cd <- cbind(cd,counts_t);
    print(i)
}

cd$test <- NULL
cd_s_2_1 <- as(as.matrix(cd), "dgCMatrix")
cd <- NULL
print("PBSFE18")

########Load files
folder_names <- c("PBS2FE18")
if (species == "h") cd <- data.frame(test=c(1:25463)) else cd <- data.frame(test=c(1:25289));

dim(cd)

for (i in folder_names)
{
    counts <- read.delim(paste(mapping_folder,i,"/",i,".counts.tsv",sep=""),skip=0,header=T,sep="\t",stringsAsFactors=F,row.names=1,check.names = F);
    rownames(counts) <- paste(i,rownames(counts),sep = "_")
    counts <- counts[rowSums(counts)>min_genes,];
    counts_t <- as.data.frame(t(counts))
    cd <- cbind(cd,counts_t);
    print(i)
}

cd$test <- NULL
cd_s_2_2 <- as(as.matrix(cd), "dgCMatrix")
cd <- NULL
print("PBS2FE18")

############Load files
folder_names <- c("MIAME18")
if (species == "h") cd <- data.frame(test=c(1:25463)) else cd <- data.frame(test=c(1:25289));

dim(cd)

for (i in folder_names)
{
    counts <- read.delim(paste(mapping_folder,i,"/",i,".counts.tsv",sep=""),skip=0,header=T,sep="\t",stringsAsFactors=F,row.names=1,check.names = F);
    rownames(counts) <- paste(i,rownames(counts),sep = "_")
    counts <- counts[rowSums(counts)>min_genes,];
    counts_t <- as.data.frame(t(counts))
    cd <- cbind(cd,counts_t);
    print(i)
}

cd$test <- NULL
cd_s_3_1 <- as(as.matrix(cd), "dgCMatrix")
cd <- NULL
print("MIAME18")

################Load files
folder_names <- c("MIA2ME18")
if (species == "h") cd <- data.frame(test=c(1:25463)) else cd <- data.frame(test=c(1:25289));

dim(cd)

for (i in folder_names)
{
    counts <- read.delim(paste(mapping_folder,i,"/",i,".counts.tsv",sep=""),skip=0,header=T,sep="\t",stringsAsFactors=F,row.names=1,check.names = F);
    rownames(counts) <- paste(i,rownames(counts),sep = "_")
    counts <- counts[rowSums(counts)>min_genes,];
    counts_t <- as.data.frame(t(counts))
    cd <- cbind(cd,counts_t);
    print(i)
}

cd$test <- NULL
cd_s_3_2 <- as(as.matrix(cd), "dgCMatrix")
cd <- NULL
print("MIA2ME18")

########Load files
folder_names <- c("MIAFE18")
if (species == "h") cd <- data.frame(test=c(1:25463)) else cd <- data.frame(test=c(1:25289));

dim(cd)

for (i in folder_names)
{
    counts <- read.delim(paste(mapping_folder,i,"/",i,".counts.tsv",sep=""),skip=0,header=T,sep="\t",stringsAsFactors=F,row.names=1,check.names = F);
    rownames(counts) <- paste(i,rownames(counts),sep = "_")
    counts <- counts[rowSums(counts)>min_genes,];
    counts_t <- as.data.frame(t(counts))
    cd <- cbind(cd,counts_t);
    print(i)
}

cd$test <- NULL
cd_s_4_1 <- as(as.matrix(cd), "dgCMatrix")
cd <- NULL
print("MIAFE18")

########Load files
folder_names <- c("MIA2FE18")
if (species == "h") cd <- data.frame(test=c(1:25463)) else cd <- data.frame(test=c(1:25289));

dim(cd)

for (i in folder_names)
{
    counts <- read.delim(paste(mapping_folder,i,"/",i,".counts.tsv",sep=""),skip=0,header=T,sep="\t",stringsAsFactors=F,row.names=1,check.names = F);
    rownames(counts) <- paste(i,rownames(counts),sep = "_")
    counts <- counts[rowSums(counts)>min_genes,];
    counts_t <- as.data.frame(t(counts))
    cd <- cbind(cd,counts_t);
    print(i)
}

cd$test <- NULL
cd_s_4_2 <- as(as.matrix(cd), "dgCMatrix")
cd <- NULL
print("MIA2FE18")


cd_comb = cbind(cd_s_1_1, cd_s_1_2, cd_s_2_1, cd_s_2_2, cd_s_3_1, cd_s_3_2, cd_s_4_1, cd_s_4_2)
dim(cd_comb)

save(cd_comb,file = "cd_s_comb.Robj")

print("file loading finished")

### Seurat
seurat_mat <- CreateSeuratObject(counts = cd_comb, min.cells = 3, project = projectName)

seurat_mat[["percent.mt"]] <- PercentageFeatureSet(seurat_mat, pattern = "^MT-")

seurat_mat[["percent.ribo"]] <- PercentageFeatureSet(seurat_mat, pattern = "^RP[SL]")

head(seurat_mat@meta.data, 5)

seurat_mat <- subset(seurat_mat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#Normalize
seurat_mat <- NormalizeData(object = seurat_mat, normalization.method = "LogNormalize", scale.factor = 10000)
print("Normalization is completed")

#ID variable genes
seurat_mat <- FindVariableFeatures(seurat_mat, selection.method = "vst", nfeatures = 2000)
print("FindVariableGenes is completed")

#Regression and Scaling with 8 cores!
seurat.mat <- ScaleData(pbmc, vars.to.regress = "percent.mt", "percent.ribo")
print("ScaleData and Regression is completed")

#Run PCA
seurat_mat <- RunPCA(seurat_mat, features = VariableFeatures(object = seurat_mat))

pdf("Elbowplot.pdf")
ElbowPlot(seurat_mat)
dev.off()

##Find clusters, VlnPlot per cluster, and see whether some of your clusters are formed of cells that systematically have lower or higher number of expressed genes
seurat_mat <- FindNeighbors(seurat_mat, dims = 1:30)
seurat_mat <- FindClusters(seurat_mat, resolution = res)

save(seurat_mat, file = "seurat_mat.Robj")

saveRDS(seurat_mat, file = "seurat_mat.rds")

#Run UMAP and save coordinates
seurat_mat <- RunUMAP(seurat_mat, dims = 1:30)

embeds = Embeddings(seurat_mat[["umap"]])

write.csv(embeds, file = paste("umapcoordinates",30,"_",res,".csv",sep=""))

seuratClusters <- seurat_mat[[c("seurat_clusters")]]

write.csv(seuratClusters, file = "seuratClusters.csv")

pdf("umap.pdf")
DimPlot(seurat_mat, reduction = "umap", labels = TRUE, pt.size = 0.5) + NoLegend()
dev.off()
print("RunUmap is done")


ClusterMarkers <- FindAllMarkers(seurat_mat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ClusterMarkers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(ClusterMarkers, file = "ClusterMarkers_subset.csv")

save(seurat_mat, file = "seurat_mat.Robj")

print("Seurat_mat has been saved")






