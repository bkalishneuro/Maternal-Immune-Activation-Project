#### Seurat Main Cluster analysis on O2
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)

#Inputs (in this case, specific to E18)
date="20200207"
species="m"
projectName ="MIA_E18"
base_dir <- "/n/scratch2/bf78/MIASingleCell/MIA18"
data_dir <- "/n/scratch2/bf78/MIASingleCell/MIA18"
setwd = "/n/scratch2/bf78/MIASingleCell/MIA18"


#Only include cells in the final seurat object
seuratClusters <- read.csv("seuratClusters_NewFinalE18.csv",header = T, row.names = 1,check.names = F)
colnames(seuratClusters)="cluster"

#load in E18 combined counts object
load(paste(data_dir,"/","cd_s_comb.Robj",sep = ""))

cd <- (cd_comb)[,rownames(seuratClusters)]

df1 <- cd[1:5000, ]
df2 <- cd[5001:10000, ]
df2 <- cd[10001:15000, ]
df2 <- cd[5001:10000, ]
df3 <- cd[10001:15000, ]
df4 <- cd[15001:20000, ]
df5 <- cd[20001:25289, ]

coverage <- colSums(cd)/10000

df1_norm <- t(apply(df1,1,function(x) x/coverage))
df2_norm <- t(apply(df2,1,function(x) x/coverage))
df3_norm <- t(apply(df3,1,function(x) x/coverage))
df4_norm <- t(apply(df4,1,function(x) x/coverage))
df5_norm <- t(apply(df5,1,function(x) x/coverage))

dfinal1 <- rbind(df1_norm, df2_norm)
dfinal2 <- rbind(dfinal1, df3_norm)
dfinal3 <- rbind(dfinal2, df4_norm)
dfinal4 <- rbind(dfinal3, df5_norm)

cd_norm <- dfinal4

#Check that all columns sum to 10000
colSums(cd_norm)

#Create normalized matrix
counts.norm_t <- data.frame(t(cd_norm),check.names = F)
counts.norm_t$seuratCluster <- factor(unlist(seuratClusters$cluster))
#Add folder and samples names appropriate for time point and experimental group
folder_names <- c("PBSFE18", "PBS2FE18", "PBSME18", "PBS2ME18", "MIAFE18", "MIA2FE18", "MIAME18", "MIA2ME18")
sample_assign <- c("PBSFE18", "PBS2FE18", "PBSME18", "PBS2ME18", "MIAFE18", "MIA2FE18", "MIAME18", "MIA2ME18")
stim_assign <- c("PBSF", "PBSF", "PBSM", "PBSM", "MIAF", "MIAF", "MIAM","MIAM")
type_assign <- c("female", "female", "male", "male", "female", "female", "male", "male")

sample_index <- data.frame(row.names = folder_names, sample=sample_assign, stim=stim_assign, type=type_assign, stringsAsFactors = F)

sample_vector <- vector(length = length(colnames(cd)));
stim_vector <- vector(length = length(colnames(cd)));
type_vector <- vector(length = length(colnames(cd)));
#replica_vector <- vector(length = length(colnames(cd)));

#This part is different in 10x and indrop. Here I say finish the cell name with -Number
for (i in rownames(sample_index))
{
  stim_vector[grep(i,colnames(cd),value = FALSE)]=sample_index[i,"stim"];
  sample_vector[grep(i,colnames(cd),value = FALSE)]=sample_index[i,"sample"];
  type_vector[grep(i,colnames(cd),value = FALSE)]=sample_index[i,"type"];
}

sample_factor<-factor(sample_vector)
stim_factor<-factor(stim_vector)
type_factor<-factor(type_vector)


### Create normalized counts table

counts.norm_t <- as.data.frame(t(cd_norm))
counts.norm_t$stim<-stim_factor
counts.norm_t$sample<-sample_factor
counts.norm_t$type<-type_factor
counts.norm_t$cluster <- seuratClusters$cluster
counts.norm_t$clusterName <- seuratClusters$clusterName

##############load UMAP####################
umap_df <- read.delim("umapcoordinates_NewFinalE18_30_1.csv",header = T,row.names=1,sep =",")
umap_df <- umap_df[rownames(seuratClusters),]
head(umap_df)
counts.norm_t <- cbind(counts.norm_t,umap_df)
seuratClusters <- cbind(seuratClusters,umap_df)
seuratClusters$sample<-sample_factor
seuratClusters$stim<-stim_factor
seuratClusters$type<-type_factor
counts.norm_t$UMI <- colSums(cd)

print("umap is loaded")

save(counts.norm_t, file = "counts.norm_t_ClusterNames_NewFinalE18.Robj")
save(stim_factor, file = "stim_factor_ClusterNames_NewFinalE18.Robj")
save(sample_factor, file = "sample_factor_ClusterNames_NewFinalE18.Robj")
save(type_factor, file = "type_factor_ClusterNames_NewFinalE18.Robj")

