https://github.com/dmcable/spacexr
https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html


# install.packages("devtools")
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)


library(spacexr)
library(Matrix)
library(data.table)
library(Seurat)

#####################################################################################
######## scRNAseq data
#####################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/BroadInstitute_SingleCell")
count = data.frame(fread("DATA_MATRIX_LOG_TPM.txt"))
rownames(count) = as.character(count$GENE)
colnames(count) = gsub("[.]","_",colnames(count))
count = count[,-1]

annot = read.delim("CLUSTER_AND_SUBCLUSTER_INDEX.txt",sep='\t',skip=1,header=T)
annot$TYPE = gsub("-","_",as.character(annot$TYPE))
annot$TYPE = gsub("[.]","_",as.character(annot$TYPE))
colnames(annot)[2] = "subclass"
rownames(annot) = as.character(annot$TYPE)

all( colnames(count) %in% as.character(annot$TYPE) )
annot = annot[ match(colnames(count), as.character(annot$TYPE) ) , ]
all( colnames(count) == rownames(annot) )
index = which( !as.character(annot$subclass)%in%c("Ependymal","Non") )

#####################################################################################
if(log) {
    counts = (count[,index])
    meta_data = annot[index,]
    resultsdir <- "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/RCTD"
} else {
    counts = 2^(count[,index])
    meta_data = annot[index,]
    resultsdir <- "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/RCTD_unlog"
}
#####################################################################################

cell_types <- meta_data$subclass
names(cell_types) <- rownames(meta_data) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- colSums(counts)
names(nUMI) <- rownames(meta_data) # create nUMI named list
#####################################################################################
reference <- Reference(counts, cell_types, nUMI, require_int = F)
print(dim(reference@counts)) #observe Digital Gene Expression matrix
table(reference@cell_types) #number of occurences for each cell type


#####################################################################################
######## spatial data
#####################################################################################
datadir <- '/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/Visium_FFPE_Mouse_Brain'
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


