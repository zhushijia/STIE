https://github.com/dmcable/spacexr
https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html


# install.packages("devtools")
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)


library(spacexr)
library(Matrix)
library(Seurat)

#####################################################################################
######## scRNAseq data
#####################################################################################

refdir <- "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq"
setwd(refdir)
counts <- readMM("count_matrix_sparse.mtx")
genes <- read.delim("count_matrix_genes.tsv",header=F)
cell_ids <- read.delim("count_matrix_barcodes.tsv",header=F)
rownames(counts) <- as.character(genes[,1])
colnames(counts) <- as.character(cell_ids[,1])

meta_data = read.csv("metadata.csv",header=T,row.names=1)
cell_types <- meta_data$celltype_major
names(cell_types) <- rownames(meta_data) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- meta_data$nCount_RNA
names(nUMI) <- rownames(meta_data) # create nUMI named list

#####################################################################################
######### downsampling to 15%
#####################################################################################
uni_celltypes = unique(cell_types)
set.seed(1234567)
tmp = lapply(uni_celltypes, function(x) {
    index = which(cell_types==x)
    sample( index, ceiling(length(index)*0.15) )
} )
selected = c()
for(k in 1:length(tmp)) selected = c(selected, tmp[[k]])
counts = counts[,selected]
cell_types = cell_types[selected]
#####################################################################################
reference <- Reference(counts, cell_types, nUMI)
print(dim(reference@counts)) #observe Digital Gene Expression matrix
table(reference@cell_types) #number of occurences for each cell type


#####################################################################################
######## spatial data
#####################################################################################
datadir <- '/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/Visium_FFPE_Human_Breast_Cancer'
matrix_paths <- paste( datadir, "outs/filtered_feature_bc_matrix.h5", sep="/")
tissue_paths <- paste( datadir, "outs/spatial/tissue_positions_list.csv", sep="/")
counts <- as.matrix(Read10X_h5(matrix_paths))
coords <- read.csv(tissue_paths, header=F, row.names=1)
coords <- coords[,c(4:5)]
colnames(coords) = c("xcoord", "ycoord")
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI

### Create SpatialRNA object
puck <- SpatialRNA(coords, counts, nUMI)

## Examine SpatialRNA object (optional)
print(dim(puck@counts)) # observe Digital Gene Expression matrix
hist(log(puck@nUMI,2)) # histogram of log_2 nUMI


print(head(puck@coords)) # start of coordinate data.frame
barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 
# This list can be restricted if you want to crop the puck e.g. 
# puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
# on the plot:
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                     title ='plot of nUMI') 


#####################################################################################
######## run RCTD
#####################################################################################

myRCTD <- create.RCTD(puck, reference, max_cores = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

results <- myRCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA

resultsdir <- "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/RCTD"
dir.create(resultsdir)
setwd(resultsdir)
save( myRCTD, norm_weights, file="myRCTD.RData" )



plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
# Plots all weights for each cell type as in full_mode. (saved as 
# 'results/cell_type_weights.pdf')
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
# Plots the weights for each cell type as in doublet_mode. (saved as 
# 'results/cell_type_weights_doublets.pdf')
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
                     results$results_df) 
# Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 
# 'results/cell_type_occur.pdf')
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)


plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir) 

# doublets
#obtain a dataframe of only doublets
doublets <- results$results_df[results$results_df$spot_class == "doublet_certain",] 
# Plots all doublets in space (saved as 
# 'results/all_doublets.pdf')
plot_doublets(spatialRNA, doublets, resultsdir, cell_type_names) 


plot_doublets_type(spatialRNA, doublets, resultsdir, cell_type_names) 
# a table of frequency of doublet pairs 
doub_occur <- table(doublets$second_type, doublets$first_type) 
# Plots a stacked bar plot of doublet ocurrences (saved as 
# 'results/doublet_stacked_bar.pdf')

plot_doub_occur_stack(doub_occur, resultsdir, cell_type_names) 


