#' calculate_spot_radius
#'
#' calculate_spot_radius calculates the spot radius based on the ratio between spot radius and spot center distance 
#'
#' @param spot_coordinates a data frame representing the spot coordinates, 
#' where the "barcode" column represents the spot id, "pixel_x" and "pixel_y" represent the coordinates on x-axis and y-axis, respectively
#' @param fct a numeric value representing the ratio between spot radius and spot center distance
#'
#' @return a numeric value representing the spot radius
#' @export
#' 
#' @author Shijia Zhu, \email{shijia.zhu@@utsouthwestern.edu}
#'
#' @seealso \code{\link{get_cells_on_spot}}; 
#' 
#' 
calculate_spot_radius <- function(spot_coordinates, fct)
{
    distMat <- dist(spot_coordinates[, c("pixel_x", "pixel_y")])
    distMat <- as.matrix(distMat)
    diag(distMat) <- Inf
    center_to_center_distance <- mean(apply(distMat, 2, min))
    
    cat("Defining spot radius...\n")
    center_to_center_distance*fct
    
}
