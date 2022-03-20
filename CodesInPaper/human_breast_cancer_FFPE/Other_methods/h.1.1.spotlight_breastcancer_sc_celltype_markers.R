

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
library(Matrix)
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq")
annot = read.csv("metadata.csv",header=T,row.names=1)
counts <- readMM("count_matrix_sparse.mtx")
genes <- read.delim("count_matrix_genes.tsv",header=F)
cell_ids <- read.delim("count_matrix_barcodes.tsv",header=F)
rownames(counts) <- as.character(genes[,1])
colnames(counts) <- as.character(cell_ids[,1])

dataSC = counts[,match( rownames(annot), colnames(counts) )]
all(rownames(annot)==colnames(dataSC))
annot$celltypes = gsub(" |-", "", as.character(annot$celltype_major) )

uni_celltypes = unique(annot$celltypes)
uni_celltypes = setdiff( uni_celltypes, c("PVL","Myeloid") )
set.seed(1234567)
selected = do.call(c, lapply(uni_celltypes, function(x) {
    index = which(celltypes==x)
    sample( index, ceiling(length(index)*0.15) )
} ))

breastcancer_sc <- CreateSeuratObject(counts=dataSC[,selected], meta.data=annot[selected,]  ) 

####################################################################################

set.seed(123)
breastcancer_sc <- Seurat::SCTransform(breastcancer_sc, verbose = FALSE) %>%
    Seurat::RunPCA(., verbose = FALSE) %>%
    Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)

Seurat::Idents(object = breastcancer_sc) <- breastcancer_sc@meta.data$celltypes
cluster_markers_all <- Seurat::FindAllMarkers(object = breastcancer_sc, 
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE)


setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/SPOTlight")
save( breastcancer_sc, cluster_markers_all, file ="breastcancer_sc_cluster_markers_all.RData"  )



