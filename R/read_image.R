
read_image= function(image, croparea=NULL, x_scale=1 ) {
    
    cat("Reading", image, "...\n")
    
    im <- image_read(path = image)
    iminf <- image_info(im)

    if (!is.null(croparea)) {
        
        cat("Cropping image to a rectangle specified by:", croparea, "...\n")
        # read crop rectangle
        rect <- as.numeric(unlist(strsplit(croparea, ",")))
        stopifnot(length(rect) == 4)
        w <- rect[1]; h <- rect[2]; xoff <- rect[3]; yoff <- rect[4]
        
        im <- im %>%
            image_crop(geometry = geometry_area(width = w, height = h, x_off = xoff, y_off = yoff))
    }
    
    stopifnot(x_scale <= 1 & x_scale > 0)
    
    if (x_scale<1) {
        
        iminf <- image_info(im)
        print(paste0("Image width: ", iminf$width))
        print(paste0("Image width after scaling: ", x_scale*iminf$width))
        
        im <- im %>%
            image_scale(paste0(round(x_scale*iminf$width))) 
        
    }

    im <- im %>%
        as_EBImage()
    
    im
}
