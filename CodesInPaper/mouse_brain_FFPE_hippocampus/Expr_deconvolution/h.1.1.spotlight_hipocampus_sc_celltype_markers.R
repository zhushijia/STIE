

module load R/4.1.1-gccmkl
module load hdf5

library(Matrix)
library(data.table)
library(Seurat)
library(SeuratData)
library(dplyr)
library(gt)
library(igraph)
library(RColorBrewer)
library(hdf5r)
library(SPOTlight)

####################################################################################
########### single cell
####################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/BroadInstitute_SingleCell")
count = data.frame(fread("DATA_MATRIX_LOG_TPM.txt"))
rownames(count) = as.character(count$GENE)
colnames(count) = gsub("[.]","_",colnames(count))
count = count[,-1]
#count = exp(count) - 1
annot = read.delim("CLUSTER_AND_SUBCLUSTER_INDEX.txt",sep='\t',skip=1,header=T)
annot$TYPE = gsub("-","_",as.character(annot$TYPE))
annot$TYPE = gsub("[.]","_",as.character(annot$TYPE))
colnames(annot)[2] = "subclass"
rownames(annot) = as.character(annot$TYPE)

all( colnames(count) %in% as.character(annot$TYPE) )
annot = annot[ match(colnames(count), as.character(annot$TYPE) ) , ]
all( colnames(count) == rownames(annot) )

index = which( !as.character(annot$subclass)%in%c("Ependymal","Non") )

hipo_sc <- CreateSeuratObject(counts=count[,index], meta.data=annot[index,] ) 
####################################################################################

set.seed(123)
hipo_sc <- Seurat::SCTransform(hipo_sc, verbose = FALSE) %>%
    Seurat::RunPCA(., verbose = FALSE) %>%
    Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)

Seurat::Idents(object = hipo_sc) <- hipo_sc@meta.data$subclass
cluster_markers_all <- Seurat::FindAllMarkers(object = hipo_sc, 
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE)


outputDir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/SPOTlight"
dir.create(outputDir, recursive = T)
setwd(outputDir)
save( hipo_sc, cluster_markers_all, file ="hipo_sc_cluster_markers_all.RData"  )



