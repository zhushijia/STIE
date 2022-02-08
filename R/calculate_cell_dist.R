#' calculate_cell_dist
#'
#' @param cells_on_spot a data frame indicating the cells on spots
#' @param max_dist a numeric value indicating the cutoff of maximum distance between two cells
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
#' 
#' @author Shijia Zhu, \email{shijia.zhu@@UTsouthwestern.edu}
#'
#' @references
#'
#' @seealso \code{\link{get_cells_on_spot}}; 
#'
#'
calculate_cell_dist <- function(cells_on_spot, max_dist=Inf)
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
