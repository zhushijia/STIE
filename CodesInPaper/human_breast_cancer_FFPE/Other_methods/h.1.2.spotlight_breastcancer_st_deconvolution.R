module load R/4.1.1-gccmkl
module load hdf5

library(Matrix)
library(data.table)
library(Seurat)
library(SeuratData)
library(dplyr)
library(gt)
library(SPOTlight)
library(igraph)
library(RColorBrewer)
library(hdf5r)

####################################################################################
############## load single cell cell types markers
####################################################################################

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/SPOTlight")
load("breastcancer_sc_cluster_markers_all.RData")

####################################################################################
############## spatial transcriptome
####################################################################################
data_dir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/Visium_FFPE_Human_Breast_Cancer/outs"
breastcancer_st = Load10X_Spatial(data.dir = data_dir)

#0.6.1 SPOTlight Decomposition
set.seed(123)

spotlight_ls <- spotlight_deconvolution(
    se_sc = breastcancer_sc,
    counts_spatial = breastcancer_st@assays$Spatial@counts,
    clust_vr = "celltypes", # Variable in sc_seu containing the cell-type annotation
    cluster_markers = cluster_markers_all, # Dataframe with the marker genes
    cl_n = 100, # number of cells per cell type to use
    hvg = 3000, # Number of HVG to use
    ntop = NULL, # How many of the marker genes to use (by default all)
    transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
    method = "nsNMF", # Factorization method
    min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
)

saveRDS(object = spotlight_ls, file = "spotlight_ls.rds")

#Read RDS object
#spotlight_ls <- readRDS(file = "spotlight_ls.rds" )
nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]

#0.7 Visualization
decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
rownames(decon_mtrx) <- colnames(breastcancer_st)

decon_df <- decon_mtrx %>%
    data.frame() %>%
    tibble::rownames_to_column("barcodes")

breastcancer_st@meta.data <- breastcancer_st@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(decon_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")

cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]
