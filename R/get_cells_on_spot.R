get_cells_on_spot <- function( cell_coordinates, spot_coordinates, spot_radius)
{
    set1 <- spot_coordinates[, c("pixel_x", "pixel_y")]
    set2 <- cell_coordinates[, c("pixel_x", "pixel_y")]
    mindist <- apply(set1, 1, function(x) {
        sqrt(colSums((t(set2[, 1:2]) - x)^2))
    })
    #spotids <- paste0(spot_coordinates$x, "x", spot_coordinates$y)
    spotids <- as.character(spot_coordinates$barcode)
    
    arrind <- which(mindist <= spot_radius, arr.ind = TRUE)
    df <- data.frame(dist = apply(arrind, 1, function(x) {mindist[x[1], x[2]]}),
                     ID = arrind[, "row"], spot = spotids[arrind[, "col"]])
    cell_coordinates$within <- "no"
    cell_coordinates$within[df$ID] <- "yes"
    cell_coordinates$spot <- NA
    cell_coordinates$spot[df$ID] <- as.character(df$spot)
    
    cat("Summarizing cells per spot ...\n")
    
    cells_on_spot <- na.omit(cell_coordinates)
    #cells_on_spot$Eccentricity <- with(cells_on_spot, sqrt( 1-Minor^2/Major^2 ) )
    
    cells_on_spot
}

