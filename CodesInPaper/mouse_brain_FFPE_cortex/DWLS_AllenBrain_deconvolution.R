
library(Seurat)
library(dplyr)
library(hdf5r)
library(data.table)
source("/archive/SCCC/Hoshida_lab/s184554/Project/stRNAseq/DWLS/Deconvolution_functions.R")

workdir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/allen_cortex"
setwd(workdir)
x = load("AllenCortex_MouseBrain_scRNASeq_DWLS_Signature.RData")

matrix_paths = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/Visium_FFPE_Mouse_Brain/outs/filtered_feature_bc_matrix.h5"
dataST <- as.matrix(t(Read10X_h5(matrix_paths)))

Signature = Signature[rownames(Signature)%in%colnames(dataST),]
Signature = Signature/100


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
save(prop_mat, file="Visium_FFPE_Mouse_Brain_cortex_DWLS_on_allen_cortex_Signature.RData" )
















