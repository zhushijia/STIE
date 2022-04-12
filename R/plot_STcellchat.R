#' plot_STcellchat
#'
#' @param cellchat 
#' @param cell_dist 
#' @param im 
#' @param image_transparency 
#' @param w 
#' @param h 
#' @param xoff 
#' @param yoff 
#' @param x_scale 
#' @param plot_object 
#' @param direction 
#' @param color_use 
#' @param lwd 
#' @param plot_cell 
#' @param contour 
#'
#' @return
#' @export
#'
#' @examples
plot_STcellchat <- function(cellchat, cell_dist, 
                            im, image_transparency=0, 
                            w=NULL, h=NULL, xoff=0, yoff=0, x_scale=1,
                            plot_object, direction=c("outgoing","incoming"), 
                            color_use=NULL, lwd=10,
                            plot_cell=FALSE, contour=NULL) {
    
    # 1. if pathway, draw on pathway
    # 2. else draw weight of all pathways
    # 3. draw pattern together and separately
    # 4. draw in and out
    
    direction = match.arg(direction)
    
    ######### pathway_probs
    # count <- cellchat@net$count
    # weight <- cellchat@net$weight
    pathway_names = cellchat@netP$pathways
    pathway_probs = cellchat@netP$prob
    
    patternSignaling <- methods::slot(cellchat, "netP")$pattern [[ direction ]]
    data = patternSignaling$data
    data1 = patternSignaling$pattern$cell
    data2 = patternSignaling$pattern$signaling
    
    ######### pattern_probs
    pattern_names = sort(unique(as.character(data1$Pattern)))
    pattern_probs = lapply( pattern_names, function(pattern) {
        tmp = lapply( pathway_names, function(pathway) {
            r2 = subset(data2, Pattern==pattern & Signaling==pathway)$Contribution    
            pathway_probs[,,pathway] * r2
        })
        dat = tmp[[1]]
        for(i in 2:length(tmp)) dat=dat+tmp[[i]]
        dat
    })
    names(pattern_probs) = pattern_names
    
    if( ! plot_object%in%pattern_names & ! plot_object%in%pathway_names )
    {
        cat( "ERROR: plot_object should be in the followings:\n" )
        cat( "Patterns:\n" )
        print(pattern_names)
        cat( "Pathways:\n" )
        print(pathway_names)
        
    } else {
        
        if( !is.null(w) & !is.null(h) ) {
            im <- im %>% image_crop(geometry = geometry_area(width = w, height = h, x_off = xoff, y_off = yoff))
        } else {
            iminf <- image_info(im)
            w = iminf$width
            h = iminf$height
        }
        
        if (x_scale<1) {
            iminf <- image_info(im)
            print(paste0("Image width: ", iminf$width))
            print(paste0("Image width after scaling: ", x_scale*iminf$width))
            im <- im %>%
                image_scale(paste0(round(x_scale*iminf$width))) 
        }
        
        im <- im %>% as_EBImage()
        
        if( image_transparency>=0 & image_transparency<=1 ) 
        {
            im <- (1-image_transparency)*im + image_transparency*(im*0+1)
        }
        
        ########################################################################################
        cells_on_spot = STIE_result$cells_on_spot
        cell_types = STIE_result$cell_types
        
        uni_celltypes = sort(unique(cell_types))
        if(is.null(color_use)) color_use = get_my_colors(length(uni_celltypes), mode=2)
        names(color_use) = uni_celltypes
        cell_cols = color_use[cell_types]
        #interaction_cols = color_use[cell_types]
        
        pixel_x = ( cells_on_spot$X - xoff )*x_scale
        pixel_y = ( cells_on_spot$Y - yoff )*x_scale
        cell_coordinates = data.frame(pixel_x, pixel_y)
        
        #################################################################################
        ######### plot interaction
        #################################################################################
        
        if(plot_object%in%pattern_names) cci = pattern_probs[[ plot_object ]]
        if(plot_object%in%pathway_names) cci = pathway_probs[ , , plot_object ]
        
        mfrow = par()$mfrow
        mar = par()$mar
        
        mat1 <- matrix(1:2, ncol=2, nrow=1)
        layout(mat1, widths=c(8, 2), heights=c(1, 2))
        
        plot(im)
        
        cat("Plotting interactions ... \n")
        for(k in 1:nrow(cell_dist))
        {
            #cat(k,"\n")
            i = cell_dist$i[k]
            j = cell_dist$j[k]
            ci = cell_types[i]
            cj = cell_types[j]
            if(direction=="outgoing") lines( cell_coordinates[c(i,j),], col=color_use[ci], lwd=lwd*cci[ci,cj] )
            if(direction=="incoming") lines( cell_coordinates[c(i,j),], col=color_use[cj], lwd=lwd*cci[cj,ci]  )
        }
        
        if(plot_cell) {
            contour = contour[ match(names(cell_types), names(contour)) ]
            plot_cell_contour(contour, cells_on_spot, w, h, xoff, yoff, x_scale, cell_cols)
        }
        
        par(mar = c(0, 0, 0, 0))
        plot.new()
        legend = uni_celltypes
        legend('topleft', legend=legend, lty=1, lwd=3, col=color_use, 
               box.lwd = 0, box.col = "white",bg = "white" )
        
        
        par( mfrow=mfrow )
        par( mar=mar)
        
    }
    
    
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
            if( any(x>0) & any(x<w) & any(y>0) & any(y<h) ) lines(x , y , col=cell_cols[i])#
        }
    }
    
}

