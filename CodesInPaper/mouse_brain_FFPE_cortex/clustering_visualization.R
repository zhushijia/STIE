result_batch_visualization <- function(STIE_result, im, spot_coordinates, contour, colors, outputDir,
                                       whole_scale, sub_scale=1, 
                                       w=3000, h=3000, margin=100) 
    
{
    dir.create(outputDir, recursive = TRUE)
    setwd(outputDir)
    
    PME_uni_cell = STIE_result$PME_uni_cell
    cell_types = STIE_result$cell_types
    contour2 = contour[ match(rownames(PME_uni_cell),names(contour)) ]
    
    pdf( "whole.pdf" )
    
    plot_sub_image(im=im, image_transparency=0,
                   x_scale=whole_scale, spot_coordinates=spot_coordinates, 
                   contour=contour2, cell_types=cell_types, color_use=colors, 
                   plot_spot=F, plot_cell=T, 
                   axis_tick=0, axis_col='grey'  )
    
    
    uni_celltypes = sort(unique(cell_types))
    for(j in 1:length(uni_celltypes) )
    {
        index = which(cell_types==uni_celltypes[j])
        plot_sub_image(im=im, image_transparency=0.5,
                       x_scale=whole_scale, spot_coordinates=spot_coordinates, 
                       contour=contour2[index], cell_types=cell_types[index], color_use=colors[j], 
                       plot_spot=F, plot_cell=T, 
                       axis_tick=0, axis_col='grey'  )
        
    }
    
    dev.off()
    
    iminf <- image_info(im)
    width <- iminf$width
    height <- iminf$height
    
    w_num <- ceiling(width/w)
    h_num <- ceiling(height/h)
    
    for( i in 1:w_num ) {
        
        for(j in 1:h_num) {
            
            xoff_ij <- w*(i-1)+1
            yoff_ij <- h*(j-1)+1
            
            w_ij <- w+margin
            h_ij <- h+margin
            
            cat(i, "and", j, "->  xoff:", xoff_ij, "yoff:", yoff_ij,
                "w:", w_ij, "h:", h_ij, "...\n")
            
            if( is.null(spot_coordinates) ) {
                spot_covered = TRUE
            } else {
                mm = w_ij*0.1
                spot_covered = any( mapply( function(x,y) 
                    (xoff_ij-mm) <= x & (xoff_ij+w_ij+mm) >= x & (yoff_ij-mm) <= y & (yoff_ij+h_ij+mm) >= y , 
                    as.integer(spot_coordinates$imagecol), as.integer(spot_coordinates$imagerow) ) )
            }
            
            if(spot_covered)
            {
                outpath <- paste0("split_im_w_",w_ij,"_h_",h_ij,"_xoff_",xoff_ij,"_yoff_",yoff_ij,".pdf")
                pdf(outpath)
                
                plot_sub_image(im=im, image_transparency=0,w=w_ij, h=h_ij, xoff=xoff_ij, yoff=yoff_ij, 
                               x_scale=sub_scale, spot_coordinates=spot_coordinates, 
                               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=F  )
                
                plot_sub_image(im=im, image_transparency=1,w=w_ij, h=h_ij, xoff=xoff_ij, yoff=yoff_ij, 
                               x_scale=sub_scale, spot_coordinates=spot_coordinates, 
                               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )
                

                
                
                dev.off()
            }
            
        }
    }
}


clusteringDir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/STIE"
setwd(clusteringDir)
load("MouseBrainCortex_clustering.RData")

for(i in 3:10)
{
    STIE_result = result[[i]]
    contour = cell_info$cell_contour
    outputDir = paste0( clusteringDir, "/clustering_visualization/clustering_", i )
    whole_scale=args$x_scale
    sub_scale=1
    
    myCol = c( "red", "blue", "green", "black", "cyan", "yellow", "steelblue",
               "darkorange", "darkred", "grey", "blue4", "purple", "chartreuse4", "burlywood1", "darkgoldenrod4" )
    colors = myCol[ 1:ncol(result[[i]]$Signature) ]
    cell_types = result[[i]]$cell_types
    colors = colors[ length(colors)+1-rank(table(cell_types)) ]
    
    result_batch_visualization(STIE_result=result[[i]], im=im, spot_coordinates=spot_coordinates, 
                               contour=contour, colors=colors, outputDir=outputDir,
                               whole_scale=whole_scale, sub_scale=sub_scale, w=3000, h=3000, margin=100) 
}





