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

meanExpr = t( apply(dataSC[,selected], 1, function(x) tapply(x,celltypes[selected],mean) ) )
save(meanExpr, file="Wu_etal_2021_BRCA_scRNASeq_meanExpr.RData")




####################################################################################################

library(Seurat)
library(dplyr)
library(hdf5r)
library(data.table)
source("/archive/SCCC/Hoshida_lab/s184554/Project/stRNAseq/DWLS/Deconvolution_functions.R")

workdir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/Visium_FFPE_Human_Breast_Cancer/DWLS_on_Wu_etal_2021_BRCA_Signature"
setwd(workdir)
x = load("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq/Wu_etal_2021_BRCA_scRNASeq_DWLS_Signature.RData")
Signature = get(x)
dataST <- as.matrix(t(Read10X_h5("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5")))

n_spots = nrow(dataST)
n_types = ncol(Signature)

prop_mat <- as.data.frame(matrix(0, nrow = n_spots, ncol = n_types) )
rownames(prop_mat) <- rownames(dataST)
colnames(prop_mat) <- colnames(Signature)

for (s in 1:n_spots) {
    
    print(sprintf("Estimating proportion for spot : %d / %d",
                  s,n_spots)
    )
    
    spot <- dataST[s,]
    tr <- trimData(Signature, spot)
    
    tr$sig <- tr$sig[,colSums(tr$sig) > 0]
    is_pd <- eigen(t(tr$sig)%*%tr$sig)$values
    is_pd <- all(is_pd > 10e-6)
    
    if (!(is_pd)) { 
        next
    }
    
    try(solDWLS <- solveDampenedWLS(tr$sig,tr$bulk),next)
    print("Proportions >> ")
    prop_mat[s,names(solDWLS)] <- solDWLS
}

setwd(workdir)
save(prop_mat, file="Visium_FFPE_Human_Breast_Cancer_DWLS_on_Wu_etal_2021_BRCA_Signature.RData" )


















