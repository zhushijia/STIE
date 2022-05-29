library(readbitmap)

library(grid)
library(RColorBrewer)
library(cowplot)

library(Seurat)
library(hdf5r)
library(Matrix)

library(data.table)

library("argparse")
library("data.table")
library("magick")
library("magrittr")
library("EBImage")
library("ggplot2")
library("dplyr")
library("quadprog")

library("foreach") 

library("parallel") 
library("doParallel")


library("STIE")

parentDir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HDST/SCP420"
images = list.files(parentDir, pattern="png")

image_path = paste0(parentDir, "/", images[i])
feature_dir = paste0(parentDir, "/", paste0(gsub(".png", "",images[i]),"_split_image") )
dir.create(feature_dir)


im <- image_read(image_path)

split_image(image=image_path, split_dir=feature_dir, 
            w=3000, h=3000, margin=100, x_scale=1 )

run_imageJ_plugin(
    imageJ = "/archive/SCCC/Hoshida_lab/s184554/Code/github/ImageJ/Fiji.app/ImageJ-linux64",
    plugin_macro = NULL,
    split_image_dir = feature_dir, 
    pattern="png$" )

merge_feature( feature_dir )

#################################################################
load(paste0(feature_dir,"/cell_info.RData"))
contour = cell_info$cell_contour

w = image_info(im)$width
h = image_info(im)$height
x_scale = 0.5
im2 <- im %>%
       image_scale(paste0(round(x_scale*w))) %>%
       as_EBImage()
im3 <- im2*0+1


pdf( paste0( gsub(".png", "",images[i]), "_mask.pdf") )
plot(im2)
plot(im3)
for( j in 1:length(contour) ) {
    cat(j,"\n")
    x = contour[[j]]
    lines(x[,1]*x_scale, x[,2]*x_scale, col='red') 
}  

if(i==1) {wi = 3000; hi = 3000; xoffi = 6000; yoffi = 6000; }
if(i==2) {wi = 3000; hi = 3000; xoffi = 9000; yoffi = 6000; }
if(i==3) {wi = 3000; hi = 3000; xoffi = 9000; yoffi = 3000; }
cell_types = rep('1', length(contour))
plot_sub_image(im=im, image_transparency=0,
               w = wi, h = hi, xoff = xoffi, yoff = yoffi, 
               x_scale=0.5, contour=contour, cell_types=cell_types, color_use="red", 
               plot_spot=F, plot_cell=T)
dev.off()

#################################################################


