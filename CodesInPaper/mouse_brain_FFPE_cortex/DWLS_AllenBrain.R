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




####################################################################################################

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
















