
plot_sub_image= function(im=NULL, im_path=NULL, w=NULL, h=NULL, xoff=0, yoff=0, x_scale=1, 
                         plot_spot=F, fct=0.25, spot_coordinates=NULL, spot_types=NULL,
                         plot_cell=T, contour, cell_types, axis_tick=0, axis_col='grey' ) {
    
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
    
    
    myCol = c( "green", "red", "black", "cyan", "yellow", "darkorange", "purple", "steelblue",
               "blue", "darkred", "grey", "blue4", "chartreuse4", "burlywood1", "darkgoldenrod4" )
    
    uni_celltypes = sort( unique(cell_types) )
    n = length(uni_celltypes)
    set.seed(1234)
    col = myCol[1:n]
    
    cell_cols = rep(0,length(cell_types))
    for(i in 1:n) cell_cols[cell_types==uni_celltypes[i]] = col[i]
    
    spot_cols = rep(0,length(spot_types))
    for(i in 1:n) spot_cols[spot_types==uni_celltypes[i]] = col[i]
    
    display(im,method='raster')
    
    if(plot_spot) plot_spot_info(spot_coordinates, xoff, yoff, x_scale, barcode=T, fct=fct) 
    if(plot_cell) 
    {
        plot_cell_contour(contour, cell_coordinates, w, h, xoff, yoff, x_scale, cell_cols)
        if(1) plot.new()
        legend('topright', legend=uni_celltypes, pch=1, col=col )
    }
    
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
        lines(r*cos(angles) + x, r*sin(angles) + y, col=c)
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
        plot_circle(x, y, r, 'red')
        
        #text(x, y-r, as.character(spot_coordinates$barcode[i]) )
        text(x, y-r, rownames(spot_coordinates)[i])
    }
    
}

