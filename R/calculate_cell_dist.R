#' calculate_cell_dist
#'
#' calculate_cell_dist calculates the distance between cells based on the cell coordinates
#'
#' @param cells_on_spot a data frame representing the cells on spots along with the cellular morphological features
#' @param dist_type a boolean value indicating whether the distance calculated between cell centers or cell boundary 
#' @param dist_cutoff a numeric value indicating the cutoff of maximum distance between two cells
#' @param axis a boolean value indicating whether the distance is calculated along major or minor axis
#'
#' @return a matrix of numeric values representing the distance between each cell pair
#' @export
#'
#' 
#' @author Shijia Zhu, \email{shijia.zhu@@UTsouthwestern.edu}
#'
#' @references
#'
#' @seealso \code{\link{get_cells_on_spot}}; 
#'
#'
calculate_cell_dist <- function(cells_on_spot, 
                                dist_type=c("center","boundary"), 
                                dist_cutoff=Inf,
                                axis=c("Minor","Major") )
{
    dist_type = match.arg(dist_type)
    axis = match.arg(axis)
    
    if(dist_type=="center") {
        
        X = do.call(rbind, lapply( 1:nrow(cells_on_spot), function(i) {   
            # cat(i,"\n")
            d = sqrt( ( cells_on_spot[i, "X"] - cells_on_spot[, "X"])^2 + 
                          ( cells_on_spot[i, "Y"] - cells_on_spot[, "Y"])^2 ) 
            
            j = which(d<=dist_cutoff)
            
            if(length(j)>0) j = j[ as.character(cells_on_spot$cell_id[j])!=as.character(cells_on_spot$cell_id[i]) ]
            
            res = NULL
            
            if(length(j)>0) res = data.frame(i=i,j=j,d=d[j])
            
            res
        }))
    }
    
    if(dist_type=="boundary") {
        
        X = do.call(rbind, lapply( 1:nrow(cells_on_spot), function(i) {   
            # cat(i,"\n")
            d = sqrt( ( cells_on_spot[i, "X"] - cells_on_spot[, "X"])^2 + 
                          ( cells_on_spot[i, "Y"] - cells_on_spot[, "Y"])^2 ) 
            
            dist_cutoff2 = cells_on_spot[i, axis] + cells_on_spot[, axis] + dist_cutoff 
            
            j = which( d <= dist_cutoff2   )
            
            if(length(j)>0) j = j[ as.character(cells_on_spot$cell_id[j])!=as.character(cells_on_spot$cell_id[i]) ]
            
            res = NULL
            
            if(length(j)>0) {
                dj = d[j] - cells_on_spot[j, axis] - cells_on_spot[i, axis] 
                res = data.frame(i=i,j=j,d=dj)
            }
            
            res
        }))
    }
    
    tag = apply(X,1,function(x) paste(sort(x[1:2]),collapse="_") )
    index = tapply( 1:nrow(X), tag, function(x) x[1])
    cell_dist = X[index, ]
    
    cell_dist
}


calculate_cell_dist0 <- function(cells_on_spot, max_dist=Inf)
{
    do.call(rbind, lapply( 1:nrow(cells_on_spot), function(i) {   
        cat(i,"\n")
        d = sqrt( ( cells_on_spot[i, "pixel_x"] - cells_on_spot[, "pixel_x"])^2 + 
                      ( cells_on_spot[i, "pixel_y"] - cells_on_spot[, "pixel_y"])^2 ) 
        j = which(d<=max_dist)
        j = setdiff(j,i)
        res = NULL
        if(length(j)>0) res = data.frame(i=i,j=j,d=d[j])
        res
    }))
}



