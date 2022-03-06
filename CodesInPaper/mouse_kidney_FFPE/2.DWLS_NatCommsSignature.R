library(dplyr)
library(Seurat)
library(cowplot)
library(data.table)


library(Matrix)
source("/archive/SCCC/Hoshida_lab/s184554/Project/stRNAseq/DWLS/Deconvolution_functions.R")

workdir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseKidney_FFPE/Signature/NatComms_2021"
setwd(workdir)
count = readRDS("GSE157079_P0_adult_counts.rds")
pheno = data.frame(fread("GSE157079_P0_adult_pheno.txt.gz",header=T))
cluster = data.frame(fread("GSE157079_P0_adult_clusters.txt.gz",header=T))
all( as.character(pheno$barcodes)==as.character(cluster$barcodes) )
all( as.character(pheno$barcodes)==colnames(count) )

index = which( as.character(pheno$pheno)=="Adult" & !as.character(cluster$clusters)%in%c("Unnamed","Stroma 1") )

dataSC = as.matrix(count[,index])
celltypes = as.character(cluster$clusters)[index]
celltypes = gsub(" ","",celltypes)

uni_celltypes = unique(celltypes)
set.seed(1234567)
selected = do.call(c, lapply(uni_celltypes, function(x) {
    index = which(celltypes==x)
    if(length(index)>3000) index = sample( index, 3000 )
    index
} ))

setwd(workdir)
dir.create(file.path(workdir, "DWLS_Signature"), showWarnings = FALSE)

Signature <- buildSignatureMatrixMAST(scdata=dataSC, id=celltypes, path="DWLS_Signature",
                                      diff.cutoff=0.5, pval.cutoff=0.01)
save(Signature, file="NatComms2021_AdultMouseKidney_scRNASeq_DWLS_Signature.RData")





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


















