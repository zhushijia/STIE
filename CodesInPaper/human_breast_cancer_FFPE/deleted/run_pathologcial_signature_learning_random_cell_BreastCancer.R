deconvolution = FALSE
clustering = FALSE
signature_learning = TRUE

known_signature = FALSE
known_cell_types = TRUE

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

############################################################
## mask
############################################################

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/labelling")
img = readImage("tissue_lowres_image.png")

mask_files = list.files(".",pattern="tif")
masks = lapply( mask_files, function(x) {
    y = readImage(x)
    as.matrix(y[,,1])
})
names(masks) = gsub( "Mask_|.tif", "", mask_files )
num = sapply(masks, function(x) sum(x==0) )
masks = masks[ order(num,decreasing=T) ]
masks = masks[!names(masks)%in%c("Fat","FibrousTissue")]

spot_labels = rep("unknown", nrow(spot_coordinates))
names(spot_labels) = as.character(spot_coordinates$barcode)
spot_tag = paste( round(spot_coordinates$pixel_x), round(spot_coordinates$pixel_y), sep="_" )

for( i in 1:length(masks) )
{
    cat(i,"\n")
    ind = which(masks[[i]]==0, arr.ind=T)
    ind_tag = paste(ind[,1],ind[,2],sep="_")
    spot_labels[ spot_tag%in%ind_tag ] = names(masks)[i]
}

spot_cols = rep( 'white', nrow(spot_coordinates))
spot_cols[ spot_labels=="InvasiveCarcinoma" ] = 'blue'
spot_cols[ spot_labels=="ImmuneCells" ] = 'green'
spot_cols[ spot_labels=="FibrousTissue" ] = 'purple'
spot_cols[ spot_labels=="Necrosis" ] = 'yellow'
spot_cols[ spot_labels=="Fat" ] = 'red'

plot(NA, xlim=c(100,500), ylim=c(100,500))
with( spot_coordinates, points( pixel_x , pixel_y, pch=16, col=spot_cols) )
with( spot_coordinates, points( pixel_x , pixel_y) )


############################################################
## mask
############################################################

cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2*spot_radius)
cells_on_spot$celltypes = NA
for(i in 1:length(spot_labels)) {
    cells_on_spot$celltypes[ cells_on_spot$spot==names(spot_labels)[i] ] = spot_labels[i]
}

result = STIE(ST_expr, Signature=NULL, cells_on_spot, features, lambda=0, steps=20, 
           known_signature=FALSE, known_cell_types=TRUE)



#setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
#save(res, file="HumanBreastCancer_labelling.RData")
#

cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
colors = c( "steelblue", "darkred", "black", "cyan", "yellow", "darkorange", "#4DAF4A")
colors = c( "purple", "green", "blue", "yellow")

cols[ spot_labels=="InvasiveCarcinoma" ] = 'blue'
cols[ spot_labels=="ImmuneCells" ] = 'green'
cols[ spot_labels=="FibrousTissue" ] = 'purple'
cols[ spot_labels=="Necrosis" ] = 'yellow'
cols[ spot_labels=="Fat" ] = 'red'

plot_sub_image(im=im, 
               x_scale=args$x_scale, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, 
               plot_spot=F, plot_cell=T, 
               axis_tick=0, axis_col='grey'  )


plot_sub_image(im=im, w=3000, h=3000, xoff=10000, yoff=10000, 
               x_scale=1, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

#### whole image

