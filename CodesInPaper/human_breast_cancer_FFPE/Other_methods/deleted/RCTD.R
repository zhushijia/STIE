https://github.com/dmcable/spacexr
https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html


# install.packages("devtools")
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)


library(spacexr)
library(Matrix)

refdir <- system.file("extdata",'Reference/Vignette',package = 'spacexr') # directory for the reference
counts <- read.csv(file.path(refdir,"dge.csv")) # load in counts matrix
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
meta_data <- read.csv(file.path(refdir,"meta_data.csv")) # load in meta_data (barcodes, clusters, and nUMI)
cell_types <- meta_data$cluster; names(cell_types) <- meta_data$barcode # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- meta_data$nUMI; names(nUMI) <- meta_data$barcode # create nUMI named list
reference <- Reference(counts, cell_types, nUMI)

print(dim(reference@counts)) #observe Digital Gene Expression matrix
table(reference@cell_types) #number of occurences for each cell type



datadir <- system.file("extdata",'SpatialRNA/Vignette',package = 'spacexr') # directory for sample Slide-seq dataset
counts <- read.csv(file.path(datadir,"MappedDGEForR.csv")) # load in counts matrix
coords <- read.csv(file.path(datadir,"BeadLocationsForR.csv"))
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
rownames(coords) <- coords$barcodes; coords$barcodes <- NULL # Move barcodes to rownames
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




myRCTD <- create.RCTD(puck, reference, max_cores = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')


results <- myRCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA
resultsdir <- 'RCTD_Plots' ## you may change this to a more accessible directory on your computer.
dir.create(resultsdir)
#> Warning in dir.create(resultsdir): 'RCTD_Plots' already exists


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


