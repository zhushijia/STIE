#AllenBrain Ground truth:
#http://atlas.brain-map.org/atlas?atlas=1&plate=100960252#atlas=1&plate=100960092&resolution=6.98&x=8101.333414713542&y=2466.666666666667&zoom=-2&structure=847&z=6

# Signature download
#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html#Subset_ST_for_cortex
#First dowload the seurat data from: https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1 to folder

library(dplyr)
library(Seurat)
library(cowplot)
library(data.table)


library(Matrix)
source("/archive/SCCC/Hoshida_lab/s184554/Project/stRNAseq/DWLS/Deconvolution_functions.R")

workdir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/allen_cortex"
setwd(workdir)
allen_reference <- readRDS("allen_cortex.rds")
#Idents(allen_reference) <- allen_reference$subclass

all( allen_reference@meta.data$subclass == allen_reference$subclass )

dataSC = as.matrix(allen_reference$RNA@counts)
celltypes = allen_reference@meta.data$subclass
celltypes = gsub(" ","",celltypes)
celltypes = gsub("/","or",celltypes)

setwd(workdir)
dir.create(file.path(workdir, "DWLS_Signature"), showWarnings = FALSE)

Signature <- buildSignatureMatrixMAST(scdata=dataSC, id=celltypes, path="DWLS_Signature",
                                      diff.cutoff=0.5, pval.cutoff=0.01)
save(Signature, file="AllenCortex_MouseBrain_scRNASeq_DWLS_Signature.RData")











