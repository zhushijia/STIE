#' plot_sub_image
#'
#' @param im 
#' @param im_path 
#' @param image_transparency 
#' @param w 
#' @param h 
#' @param xoff 
#' @param yoff 
#' @param x_scale 
#' @param plot_spot 
#' @param fct 
#' @param spot_coordinates 
#' @param spot_types 
#' @param plot_cell 
#' @param contour 
#' @param cell_types 
#' @param cell_types_of_interest 
#' @param color_use 
#' @param axis_tick 
#' @param axis_col 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
plot_sub_image <- function(im=NULL, im_path=NULL, image_transparency=1, w=NULL, h=NULL, xoff=0, yoff=0, x_scale=1, 
                         plot_spot=F, fct=0.25, spot_coordinates=NULL, spot_types=NULL,
                         plot_cell=T, contour, cell_types, cell_types_of_interest=unique(cell_types), 
                         color_use=NULL, axis_tick=0, axis_col='grey' ) {
    
   if( is.null(im) ) {
        cat("Reading", image, "...\n")
        im <- image_read(path = im_path)
    }
    
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

    im <- im %>%
        as_EBImage()
    
    if( image_transparency>=0 & image_transparency<=1 ) 
    {
        im <- (1-image_transparency)*im + image_transparency*(im*0+1)
    }
    
    uni_celltypes = sort( unique(cell_types) )
    num_celltypes = table(cell_types)
    num_celltypes = num_celltypes[ match(uni_celltypes, names(num_celltypes)) ]
    n = length(uni_celltypes)
    
    if(is.null(color_use))
    {
        myCol = c( "green", "red", "black", "cyan", "yellow", "darkorange", "purple", "steelblue",
                   "blue", "darkred", "grey", "blue4", "chartreuse4", "burlywood1", "darkgoldenrod4" )
        
        color_use = myCol[1:n]
    }
    
    cell_cols = rep(0,length(cell_types))
    for(i in 1:n) cell_cols[cell_types==uni_celltypes[i]] = color_use[i]
    
    spot_cols = rep(0,length(spot_types))
    for(i in 1:n) spot_cols[spot_types==uni_celltypes[i]] = color_use[i]
    
    mat1 <- matrix(1:2, ncol=2, nrow=1)
    layout(mat1, widths=c(8, 2), heights=c(1, 2))
    
    display(im,method='raster')
    
    if(plot_spot) plot_spot_info(spot_coordinates, xoff, yoff, x_scale, barcode=F, fct=fct) 
    
    if( axis_tick>0 )
    {
        w_scaled = dim(im)[1]
        h_scaled = dim(im)[2]
        x_at = ceiling(xoff/axis_tick)*axis_tick - xoff + c( 1:ceiling(w_scaled/x_scale/axis_tick) )*axis_tick
        y_at = ceiling(yoff/axis_tick)*axis_tick - yoff + c( 1:ceiling(h_scaled/x_scale/axis_tick) )*axis_tick
        x_label = ceiling(xoff/axis_tick)*axis_tick + x_at
        y_label = ceiling(yoff/axis_tick)*axis_tick + y_at
        axis(1,pos=c(0,0),at=x_at*x_scale,label=x_label)
        axis(2,pos=c(0,0),at=y_at*x_scale,label=y_label)
        
        tmp = sapply(x_at*x_scale, function(x) abline(h=x,lty=2,col=axis_col) )
        tmp = sapply(y_at*x_scale, function(x) abline(v=x,lty=2,col=axis_col) )
    }
    
    if(plot_cell) 
    {
        plot_cell_contour(contour, cell_coordinates, w, h, xoff, yoff, x_scale, cell_cols)
        par(mar = c(0, 0, 0, 0))
        plot.new()
        legend = paste0(uni_celltypes," (", num_celltypes ,")" )
        legend('topleft', legend=legend, pch=1, col=color_use, 
               box.lwd = 0, box.col = "white",bg = "white" )
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


plot_spot_info <- function(spot_coordinates, xoff, yoff, x_scale, 
                           spot_cols='green', barcode=T, fct=0.25)
{

    if(length(spot_cols)==1) spot_cols = rep( spot_cols, nrow(spot_coordinates) )
    
    plot_circle <- function(x, y, r, c) {
        angles <- seq(0, 2*pi, length.out = 360)
        lines(r*cos(angles) + x, r*sin(angles) + y, col=c, lwd=1)
    }
    
    if( all(c("X","Y") %in% colnames(spot_coordinates)) )
    {
        spot_coordinates$pixel_x = (spot_coordinates$X-xoff)*x_scale
        spot_coordinates$pixel_y = (spot_coordinates$Y-yoff)*x_scale
    }
    
    if( all(c("imagecol","imagerow") %in% colnames(spot_coordinates)) )
    {
        spot_coordinates$pixel_x = (spot_coordinates$imagecol-xoff)*x_scale
        spot_coordinates$pixel_y = (spot_coordinates$imagerow-yoff)*x_scale
    }
    
    spot_radius <- calculate_spot_radius(spot_coordinates, fct)
    
    for(i in 1:nrow(spot_coordinates)) {
        
        x = as.numeric(spot_coordinates[i, "pixel_x"])
        y = as.numeric(spot_coordinates[i, "pixel_y"])
        r = spot_radius
        c = spot_cols[i]
        plot_circle(x, y, r, 'white')
        
        #text(x, y-r, as.character(spot_coordinates$barcode[i]) )
        if(barcode) text(x, y-r, rownames(spot_coordinates)[i])
    }
    
}

