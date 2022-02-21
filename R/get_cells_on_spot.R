#' get_cells_on_spot 
#'
#' get_cells_on_spot finds the cells falling in each spot based on the cell coordinates, spot coordinates and spot radius
#' 
#' @param cell_coordinates a matix of numeric values represensting the location of each cell, where columns "pixel_x" and "pixel_y" represent the x and y coordinates, respectively
#' @param spot_coordinates a matix of numeric values represensting the location of each spot, where columns "pixel_x" and "pixel_y" represent the x and y coordinates, respectively
#' @param spot_radius a numric value, representing the spot radius
#'
#' @return A data frame containing the follow components:
#' \itemize{
#'  \item {cell_id} character values representing the unique cell ids
#'  \item {spot} character values representing the unique spot ids
#' 
#' @export
#'
#' @examples
#'
#' 
#' 
#' 
#' @author Shijia Zhu, \email{shijia.zhu@@UTsouthwestern.edu}
#'
#' @references
#'
#' @seealso \code{\link{update_morphology_parameter}}; \code{\link{get_cells_on_spot}}; 
#'
#'
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
    if(0)
    {
        cell_coordinates$within <- "no"
        cell_coordinates$within[df$ID] <- "yes"
        cell_coordinates$spot <- NA
        cell_coordinates$spot[df$ID] <- as.character(df$spot)
        
        cells_on_spot <- na.omit(cell_coordinates)
        #cells_on_spot$Eccentricity <- with(cells_on_spot, sqrt( 1-Minor^2/Major^2 ) )
        
        cat("Summarizing cells per spot ...\n")
    }
    
    cells_on_spot <- cell_coordinates[ df$ID, ]
    cells_on_spot$spot <- df$spot
    cells_on_spot$dist <- df$dist
    
    cells_on_spot
}

