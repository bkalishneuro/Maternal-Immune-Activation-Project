#Scrublet Doublet Removal and Cluster Filtering

#First, remove all cell clusters indicating poor cell health or with overrepresentation by one replicate, and create a new seurat object

#Before you start, just make sure you have a virtual python environment with scrublet installed (pip install scrublet works great). 

#The next thing to do is generate sparse raw counts matrices for each individual sample in your dataset, and save that as ".mtx" files. Once you load up an R session with your seurat object and raw counts matrix loaded, I just split my seurat object by sample, transpose rows and columns, and then save. What you want to end up with is sparse matrices where the rows are cells and the columns are genes.

List <- SplitObject(E14_UnwantedClustersRemoved,split.by="orig.ident")
Samples <- c("MIABME14","MIACME14","MIABFE14","MIACFE14","PBSBME14","PBSCME14","PBSBFE14","PBSCFE14")
for (sample in Samples){
cells= colnames(List[[sample]])
cd_comb_cells = cd_comb[,cells]
cd_comb_cells = t(cd_comb_cells)
writeMM(cd_comb_cells,file = paste0(sample,"_Test.mtx"))
}
#(and then again for each sample)

#Next, activate your python environment, and load in the .mtx files you just made:
#What we have to do is plot a histogram of both the simulated and actual doublet data distribution, and choose our own threshold. You then use this graph to find the midpoint between the peaks of the bimodal distribution of doublet scores in the simulated histogram. That value is your threshold. Cells to the right of this value will be marked as doublets. If you don't have a bimodal distribution here, then scrublet doesn't work properly, but it is okay if your "observed transcriptomes" graph doesn't look bimodal. I attached an example image of this plot. For this data, I chose a threshold of 0.3, and we got a doublet rate of around 5 percent.

scrub.plot_histogram()
plt.savefig("scrubMIABME14TestBeforeBetterFilter.png")

import scipy.io
import scrublet as scr
import matplotlib.pyplot as plt
import numpy
Samples = ["MIABME14","MIACME14","MIABFE14","MIACFE14","PBSBME14","PBSCME14","PBSBFE14","PBSCFE14"]
for sample in Samples:
	counts_matrix = scipy.io.mmread(sample+"_Test.mtx").tocsc()
	scrub = scr.Scrublet(counts_matrix)
	doublet_scores, predicted_doublets = scrub.scrub_doublets()
	if sample == "PBSCFE14":
		predicted_doublets = scrub.call_doublets(threshold=0.22)
	elif sample == "MIABFE14":
		predicted_doublets = scrub.call_doublets(threshold=0.25)
	else:
		predicted_doublets = scrub.call_doublets(threshold=0.3)
	numpy.savetxt(sample+"FinalScrublets.txt",predicted_doublets)
	scrub.plot_histogram()
	plt.savefig(sample+"TestBeforeBetterFilter.png")

#Now you can exit python, and load R back up, and load in your 0/1 lists of doublets vs no doublet cells
#and then cross-reference my 0/1 list with my cell list to eliminate doublets

library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
load("cd_s_comb.Robj")
load("E14_UnwantedClustersRemoved.Robj")
SampleList <- SplitObject(E14_UnwantedClustersRemoved,split.by = "orig.ident")
NonDoubletCells <- vector()
Samples <- c("MIABME14","MIACME14","MIABFE14","MIACFE14","PBSBME14","PBSCME14","PBSBFE14","PBSCFE14")
for (sample in Samples){
SampleScrublets <- read.delim(paste0(sample,"FinalScrublets.txt"),header = FALSE)
SampleIndex <- cbind(SampleScrublets,colnames(SampleList[[sample]]))
for (number in 1:(nrow(SampleIndex))){
if (SampleIndex[number,1] == 0){
NonDoubletCells <- c(NonDoubletCells,as.character(SampleIndex[number,2]))
}
}
}

#Assign wanted cells to new combined counts matrix object
NewComb <- cd_comb[,NonDoubletCells]

#And then create a new seurat object using my list of cells (MIAFE18Cells). I have always re-clustered after removing doublets, because I find it makes cluster identity easier to decipher.