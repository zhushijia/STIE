calculate_spot_radius <- function(spot_coordinates, fct)
{
    distMat <- dist(spot_coordinates[, c("pixel_x", "pixel_y")])
    distMat <- as.matrix(distMat)
    diag(distMat) <- Inf
    center_to_center_distance <- mean(apply(distMat, 2, min))
    
    cat("Defining spot radius...\n")
    center_to_center_distance*fct
    
}
