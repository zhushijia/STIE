
library(Seurat)
library(dplyr)
library(hdf5r)
library(data.table)
source("/archive/SCCC/Hoshida_lab/s184554/Project/stRNAseq/DWLS/Deconvolution_functions.R")

workdir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/BroadInstitute_SingleCell"

matrix_paths = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/Visium_FFPE_Mouse_Brain/outs/filtered_feature_bc_matrix.h5"
dataST <- as.matrix(t(Read10X_h5(matrix_paths)))

setwd(workdir)
Signature = read.delim("Major_cell_types_marker_genes.txt", header=T, row.names=1)
Signature = Signature[,!colnames(Signature)%in%c("Ependymal.cells")]
Signature = as.matrix(Signature)[rownames(Signature)%in%colnames(dataST),]
colnames(Signature) = gsub("Pyramidal.neurons.|.interneurons|.like.cells|Granule.cells.", "", colnames(Signature))
Signature = Signature[,order(colnames(Signature))]

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

outputDir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/DWLS"
dir.create(outputDir)
setwd(outputDir)
save(prop_mat, file="Visium_FFPE_Mouse_Brain_hippocampus_DWLS_on_BroadInstitute_Signature.RData" )





