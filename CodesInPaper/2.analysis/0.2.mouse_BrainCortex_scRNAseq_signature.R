library(dplyr)
library(Seurat)
library(cowplot)
library(data.table)
library(Matrix)

source("/archive/SCCC/Hoshida_lab/s184554/Project/stRNAseq/DWLS/Deconvolution_functions.R")

workdir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/allen_cortex"
setwd(workdir)
allen_reference <- readRDS("allen_cortex.rds")
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


