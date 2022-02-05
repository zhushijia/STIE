#' calculate_morphology_probability
#'
#' calculate_morphology_probability calculate the morphological probability by assuming the normal distribution with parameters mu and sigma
#' 
#' @param morphology_features a matrix of numeric values in size of cells*features representing the values of morphology_features for each cell, where the column names are the morphoglical features, and row names are cell ids
#' @param feature_list a vector of characters in size of morphology features representing the morphology feature names used in the analysis
#' @param mu  a vector of numeric values in size of morphology features representing the mean of the normal distribution
#' @param sigma a vector of numeric values in size of morphology features representing the standard deviation of the normal distribution
#' @param scaled boolean value representing whether the probability is scaled
#'
#' @return a matrix of numeric values in size of cells*features representing the values of morphology feature probability for each cell
#' 
#' @export
#'
#' @examples
#' # 
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
calculate_morphology_probability <- function(morphology_features, feature_list, mu, sigma, scaled=TRUE )
{
    PM = matrix( nrow=nrow(morphology_features), ncol=nrow(mu) , data=0)
    rownames(PM) = as.character(morphology_features$cell_id)
    colnames(PM) = rownames(mu)
    
    for(t in 1:nrow(mu)) {
        PM[,t] <- apply(morphology_features[,feature_list], 1, function(x) 
            prod( dnorm( as.numeric(x), mean=as.numeric(mu[t,feature_list]), sd=as.numeric(sigma[t,feature_list]) )))
    }
    
    if(scaled) PM = t( apply(PM, 1, function(x) x/sum(x) ) )
    PM
}

