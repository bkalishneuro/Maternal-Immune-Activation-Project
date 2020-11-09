library(edgeR)
library(Seurat)
library(dplyr)
library(methods)
library(Matrix)

#Load in appropriate Seurat objecct and combined counts matrix, and split object by cluster
load("Seurat_CorticalLayers.Robj")
load("cd_s_comb.Robj")
List <- SplitObject(Seurat_CorticalLayers)

for (Cluster in names(List)){
DataFrame <- cd_comb[,colnames(List[[Cluster]])]
write.tsv(DataFrame,file=paste0("E18_",Cluster,"_RawCounts.tsv"))
Group <- vector()

#Assign sample groups to cell names, adjusting "18" to "14" for E14 timepoint
for (name in colnames(DataFrame)){
if (grepl("PBSME18_",name,fixed=TRUE)){
Group <- c(Group,1)
}
else if (grepl("PBS2ME18_",name,fixed=TRUE)){
Group <- c(Group,1)
}
else if (grepl("PBSFE18_",name,fixed=TRUE)){
Group <- c(Group,2)
}
else if (grepl("PBS2FE18_",name,fixed=TRUE)){
Group <- c(Group,2)
}
else if (grepl("MIAME18_",name,fixed=TRUE)){
Group <- c(Group,3)
}
else if (grepl("MIA2ME18_",name,fixed=TRUE)){
Group <- c(Group,3)
}
else if (grepl("MIAFE18_",name,fixed=TRUE)){
Group <- c(Group,4)
}
else if (grepl("MIA2FE18_",name,fixed=TRUE)){
Group <- c(Group,4)
}
}

#Perform differential gene expression in edgeR
edgeRCounts_All <- DGEList(counts=DataFrame,group=Group)
edgeRCounts_All <- calcNormFactors(edgeRCounts_All)
design <- model.matrix(~0+group, data=edgeRCounts_All$samples)
colnames(design) <- levels(edgeRCounts_All$samples$group)
edgeRCounts_All <- estimateDisp(edgeRCounts_All, design)
fitQLF <- glmQLFit(edgeRCounts_All, design)


qlf_MIAMPBSM <- glmQLFTest(fitQLF, contrast=c(1,0,-1,0))
FDR_MIAMPBSM <- p.adjust(qlf_MIAMPBSM$table$PValue, method="BH")
qlf_MIAMPBSM_FDR <- cbind(qlf_MIAMPBSM$table,FDR_MIAMPBSM)
write.table(qlf_MIAMPBSM_FDR, file = paste0(Cluster,"_","Multifactorial_EdgeR_RNA_E14_MIAMPBSM_QValues.tsv"), sep="\t", col.names=NA)
qlf_MIAFPBSF <- glmQLFTest(fitQLF, contrast=c(0,1,0,-1))
FDR_MIAFPBSF <- p.adjust(qlf_MIAFPBSF$table$PValue, method="BH")
qlf_MIAFPBSF_FDR <- cbind(qlf_MIAFPBSF$table,FDR_MIAFPBSF)
write.table(qlf_MIAFPBSF_FDR, file = paste0(Cluster,"_","Multifactorial_EdgeR_RNA_E14_MIAFPBSF_QValues.tsv"), sep="\t", col.names=NA)
qlf_MIAMMIAF <- glmQLFTest(fitQLF, contrast=c(0,0,-1,1))
FDR_MIAMMIAF <- p.adjust(qlf_MIAMMIAF$table$PValue, method="BH")
qlf_MIAMMIAF_FDR <- cbind(qlf_MIAMMIAF$table,FDR_MIAMMIAF)
write.table(qlf_MIAMMIAF_FDR, file = paste0(Cluster,"_","Multifactorial_EdgeR_RNA_E14_MIAMMIAF_QValues.tsv"), sep="\t", col.names=NA)
qlf_InteractionTerm <- glmQLFTest(fitQLF, contrast=c(1,-1,-1,1))
FDR_InteractionTerm <- p.adjust(qlf_InteractionTerm$table$PValue, method="BH")
qlf_InteractionTerm_FDR <- cbind(qlf_InteractionTerm$table,FDR_InteractionTerm)
write.table(qlf_InteractionTerm_FDR, file = paste0(Cluster,"_","Multifactorial_EdgeR_RNA_E14_InteractionTerm_QValues.tsv"), sep="\t", col.names=NA)
qlf_PBSMPBSF <- glmQLFTest(fitQLF, contrast=c(-1,1,0,0))
FDR_PBSMPBSF <- p.adjust(qlf_PBSMPBSF$table$PValue, method="BH")
qlf_PBSMPBSF_FDR <- cbind(qlf_PBSMPBSF$table,FDR_PBSMPBSF)
write.table(qlf_PBSMPBSF_FDR, file = paste0(Cluster,"_","Multifactorial_EdgeR_RNA_E14_PBSMPBSF_QValues.tsv"), sep="\t", col.names=NA)
}