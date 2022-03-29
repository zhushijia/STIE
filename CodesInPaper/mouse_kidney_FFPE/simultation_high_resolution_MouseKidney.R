The capture area has ~5,000 gene expression spots, each spot is ~55 microns with primers that include: Illumina TruSeq Read 1 (partial read 1 sequencing primer) 16 nt Spatial Barcode (all primers in a specific spot share the same Spatial Barcode) 12 nt unique molecular identifier (UMI)

########################################################################################
######## simulate spot size
########################################################################################
source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_kidney_FFPE/parameter_MouseKidney.R")

spot_sizes = sort( c(5, 10, 20, 30, 55), decreasing=T )

STdata_sim = list()
for( i in 1:length(spot_sizes) )
{
    ss = spot_sizes[[i]]
    cat(ss, "\n")
    t = system.time( STdata_sim[[i]] <- simulate_STdata(cells_coordinates=morphology_fts, cell_types=NULL, Signature=NULL,
                                                        spot_coordinates_ref=spot_coordinates, spot_diameter_pixel_ref=spot_radius*2, spot_size_ref=55, 
                                                        spot_size_sim=ss, x_scale=args$x_scale ) )
    print(t)
    names(STdata_sim)[i] = ss
}


setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseKidney_FFPE/count/results/STIE")
save(STdata_sim, file="MouseKidney_FFPE_STdata_sim.RData" )



########################################################################################
######## plot spot covering cells
########################################################################################
source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_kidney_FFPE/parameter_MouseKidney.R")

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseKidney_FFPE/count/results/STIE")
load("MouseKidney_FFPE_STdata_sim.RData" )


########################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseKidney_FFPE/count/results/STIE")
pdf( "MouseKidney_FFPE_STdata_sim_cellcount_distribution.pdf" )



plot_circle <- function(x, y, r, c) {
    angles <- seq(0, 2*pi, length.out = 360)
    lines(r*cos(angles) + x, r*sin(angles) + y, col=c, lwd=4)
}

plot_cell_contour = function(contour, cell_coordinates, w, h, xoff, yoff, x_scale, cell_cols )
{
    
    if( is.null(contour) )
    {
        points( (cell_coordinates$X-xoff)*x_scale, (cell_coordinates$Y-yoff)*x_scale, col=cell_cols )
        
    } else {
        
        for(i in 1:length(contour))
        {
            x = ( as.numeric(contour[[i]][,1]) - xoff )*x_scale
            y = ( as.numeric(contour[[i]][,2]) - yoff )*x_scale
            if( any(x>0) & any(x<w) & any(y>0) & any(y<h) ) lines(x , y , col=cell_cols[i], lwd=3)#
        }
    }
    
}



library(vioplot)
par(mfrow=c(2,2))
cols = c('grey','darkorange','steelblue','darkgreen','darkred')
percent = lapply(STdata_sim, function(x) x$cells_on_spot$percent )
vioplot(percent, col=cols )

par(mfrow=c(2,1))
spotfreq = lapply(STdata_sim, function(x) table(x$spot_coordinates$cell_count) )
max_xlim = max( do.call(c, lapply(spotfreq, function(x) as.integer(names(x))) ) )
max_xlim = 40
max_ylim = max( do.call(c, lapply(spotfreq, function(x) log10(x) ) ) )
plot(NA, xlim=c(1,max_xlim), ylim=c(0,max_ylim), ylab="log10 spot freq", xlab="cell count" )
sapply( 1:length(spotfreq), function(i) {
    x = log10(spotfreq[[i]])
    lines( 1:length(x), x, col=cols[i])
    points( 1:length(x), x, col=cols[i], pch=16)
} )
abline(v=2, lty=2)

plot(NA, xlim=c(1,max_xlim), ylim=c(0,max_ylim), ylab="log10 spot freq", xlab="cell count" )
legend( 'topright', pch=16, col=cols, legend=rev(names(STdata_sim)) )


########################################################################################
######## plot spot example
########################################################################################


par(mfrow=c(1,2))
spot_coordinates_ref = STdata_sim[['55']]$spot_coordinates
ind = which.max( spot_coordinates_ref$cell_count )
cc = subset( STdata_sim[['55']]$cells_on_spot, spot %in% as.character(spot_coordinates_ref[ind,"barcode"]) )
overlapped_cells = as.character(cc$cell_id)

xoff = ( spot_coordinates_ref[ind,'pixel_x'] - 2*spot_radius ) / args$x_scale
yoff = ( spot_coordinates_ref[ind,'pixel_y'] - 2*spot_radius ) / args$x_scale
w = 4*spot_radius / args$x_scale
h = 4*spot_radius / args$x_scale
img <- im %>% image_crop(geometry = geometry_area(width = w, height = h, x_off = xoff, y_off = yoff))
img <- img %>% as_EBImage()
image_transparency=0.5; img <- (1-image_transparency)*img + image_transparency*(img*0+1)
display(img,method='raster')

contour1 = cell_info$cell_contour
contour2 = contour1[ names(contour1) %in% overlapped_cells ]
#plot_cell_contour(contour2, morphology_fts, w, h, xoff, yoff, x_scale=1, cell_cols=rep('green',length(contour2)) )
plot_cell_contour(contour1, morphology_fts, w, h, xoff, yoff, x_scale=1, cell_cols=rep('black',length(contour1)) )

size_draw = 5
spot_coordinates_draw = subset( STdata_sim[[as.character(size_draw)]]$spot_coordinates, 
                                imagecol>xoff & (imagecol-xoff)<w &
                                    imagerow>yoff & (imagerow-yoff)<h )
contour_draw = subset( STdata_sim[[as.character(size_draw)]]$spot_coordinates, 
                       imagecol>xoff & (imagecol-xoff)<w &
                           imagerow>yoff & (imagerow-yoff)<h )

for(i in 1:nrow(spot_coordinates_draw)) 
{
    r = size_draw*spot_radius/55/args$x_scale
    x = spot_coordinates_draw[i,'imagecol'] - xoff
    y = spot_coordinates_draw[i,'imagerow'] - yoff
    ci = 'darkred'
    plot_circle(x, y, r, ci) 
}



##############################
display(img*0+1,method='raster')
spot_sizes = sort( c(5, 10, 20, 30, 55), decreasing=T )
cols = c('grey','darkorange','steelblue','darkgreen','darkred')
for(i in 1:length(spot_sizes)) 
{
    r = spot_sizes[i]*spot_radius/55/args$x_scale
    x = w/2
    y = h/2
    ci = cols[i]
    plot_circle(x, y, r, ci) 
}



dev.off()




