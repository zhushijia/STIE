library(dplyr)
library(Seurat)
library(cowplot)
library(data.table)
library(Matrix)

source("/archive/SCCC/Hoshida_lab/s184554/Project/stRNAseq/DWLS/Deconvolution_functions.R")
workdir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq"
setwd(workdir)
annot = read.csv("metadata.csv",header=T,row.names=1)
counts <- readMM("count_matrix_sparse.mtx")
genes <- read.delim("count_matrix_genes.tsv",header=F)
cell_ids <- read.delim("count_matrix_barcodes.tsv",header=F)
rownames(counts) <- as.character(genes[,1])
colnames(counts) <- as.character(cell_ids[,1])

dataSC = counts[,match( rownames(annot), colnames(counts) )]
all(rownames(annot)==colnames(dataSC))
celltypes = gsub(" |-", "", as.character(annot$celltype_major) )

uni_celltypes = unique(celltypes)
set.seed(1234567)
selected = do.call(c, lapply(uni_celltypes, function(x) {
    index = which(celltypes==x)
    sample( index, ceiling(length(index)*0.15) )
} ))

setwd(workdir)
dir.create(file.path(workdir, "DWLS_Signature"), showWarnings = FALSE)

Signature <- buildSignatureMatrixMAST(scdata=dataSC[,selected], id=celltypes[selected], path="DWLS_Signature",
                                                 diff.cutoff=0.5, pval.cutoff=0.01)
save(Signature, file="Wu_etal_2021_BRCA_scRNASeq_DWLS_Signature.RData")













