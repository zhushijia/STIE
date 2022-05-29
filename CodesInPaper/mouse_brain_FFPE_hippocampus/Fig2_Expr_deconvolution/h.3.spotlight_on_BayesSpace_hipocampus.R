library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
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

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/SPOTlight")
load("hipo_sc_cluster_markers_all.RData")

####################################################################################
############## BayesSpace enhanced transcriptome
####################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/BayesSpace")
load("ST_enhanced.RData")

colData = data.frame(colData(ST.enhanced))
counts = logcounts(ST.enhanced)
index=apply(counts, 1, function(x) any(is.na(x)) )
counts_2 = counts[!index,]

#0.6.1 SPOTlight Decomposition
set.seed(123)

spotlight_ls <- spotlight_deconvolution(
    se_sc = hipo_sc,
    counts_spatial = counts_2,
    clust_vr = "subclass", # Variable in sc_seu containing the cell-type annotation
    cluster_markers = cluster_markers_all, # Dataframe with the marker genes
    cl_n = 100, # number of cells per cell type to use
    hvg = 3000, # Number of HVG to use
    ntop = NULL, # How many of the marker genes to use (by default all)
    transf = "raw", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
    method = "nsNMF", # Factorization method
    min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
)

saveRDS(object = spotlight_ls, file = "spotlight_ls.rds")

#Read RDS object
#spotlight_ls <- readRDS(file = "spotlight_ls.rds")

nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]
decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
rownames(decon_mtrx) <- colnames(counts_2)
cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/BayesSpace")
pdf("ST_enhanced_all_cell_types_log3.pdf")

prop_data = decon_mtrx
coordinates = colData
x_scale=1
colors = c("magenta", "blue", "green", "black", "cyan","orange")

BayesSpaceVisium_scatterpie_overlay( prop_data, coordinates, x_scale, cell_types_all, 
                       img_path=NULL, scatterpie_alpha=1, pie_scale = 0.4, 
                       colors=colors )

dev.off()




