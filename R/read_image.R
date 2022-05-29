#' read_image
#'
#' read_image loads the region of interest from an image at a predefined rescaled size
#'
#' @param image a character value indicating the path to the image
#' @param croparea a character value in the format of "width,height,xoff,yoff" for cropping the image 
#' @param x_scale a numeric value from 0 to 1 for resizing the image
#'
#' @return an EBImage object of the cropped and resized image
#' @export
#' 
#' @author Shijia Zhu, \email{shijia.zhu@@utsouthwestern.edu}
#'
read_image <- function(image, croparea=NULL, x_scale=1 ) {
    
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
