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
library(ggplot2)

####################################################################################
############## load single cell cell types markers
####################################################################################
outputDir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/SPOTlight"
setwd(outputDir)
load("hipo_sc_cluster_markers_all.RData")

####################################################################################
############## spatial transcriptome
####################################################################################
data_dir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/Visium_FFPE_Mouse_Brain/outs"
hipo_st = Load10X_Spatial(data.dir = data_dir)

#0.6.1 SPOTlight Decomposition
set.seed(123)

spotlight_ls <- spotlight_deconvolution(
    se_sc = hipo_sc,
    counts_spatial = hipo_st@assays$Spatial@counts,
    clust_vr = "subclass", # Variable in sc_seu containing the cell-type annotation
    cluster_markers = cluster_markers_all, # Dataframe with the marker genes
    cl_n = 100, # number of cells per cell type to use
    hvg = 3000, # Number of HVG to use
    ntop = NULL, # How many of the marker genes to use (by default all)
    transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
    method = "nsNMF", # Factorization method
    min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
)

saveRDS(object = spotlight_ls, file = "spotlight_ls_log.rds")

#Read RDS object
#spotlight_ls <- readRDS(file = "spotlight_ls_log.rds" )
nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]

#0.7 Visualization
decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
rownames(decon_mtrx) <- colnames(hipo_st)

decon_df <- decon_mtrx %>%
    data.frame() %>%
    tibble::rownames_to_column("barcodes")

hipo_st@meta.data <- hipo_st@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(decon_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")

cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]

pdf("all_cell_types_log.pdf")
SPOTlight::spatial_scatterpie(se_obj = hipo_st,
                              cell_types_all = cell_types_all,
                              img_path = paste0(data_dir,"/spatial/tissue_lowres_image.png"),
                              pie_scale = 0.8)
dev.off()



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### my way for backup
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
setwd(outputDir)
pdf("all_cell_types_log3.pdf",w=15,h=15)

se_obj = hipo_st
slice <- names(se_obj@images)[1]
x_scale = se_obj@images[[slice]]@scale.factors$lowres
coordinates = data.frame(se_obj@images[[slice]]@coordinates)
prop_data = se_obj@meta.data
img_path = paste0(data_dir,"/spatial/tissue_lowres_image.png")

colors = c("magenta", "blue", "green", "black", "cyan","orange")

ST_scatterpie_overlay( prop_data, coordinates, x_scale, cell_types_all, img_path=img_path, scatterpie_alpha=1, pie_scale = 0.8, colors=colors )
dev.off()


