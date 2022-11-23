#' update_morphology_parameter
#'
#' update_morphology_parameter re-etimates the parameters for the cellular morphological model in STIE
#'
#' @param PE_on_spot a matrix of cell type probabilities for each spot estimated from the spot gene expression, where the row is the spot and column is the cell type 
#' @param PM_on_cell a matrix of cell type probabilities for each cell estimated from the cellular morphological features, where the row is the cell and column is the cell type  
#' @param cells_on_spot a data frame representing the cells on spots along with the cellular morphological features, where "cell_id" represents the unique cell id and "spot" represents the uniqe spot id
#' @param features a vector of character values, representing the morphological features used in the STIE model
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
update_morphology_parameter <- function(PE_on_spot, PM_on_cell, cells_on_spot, features)
{
    spot_id = as.character(cells_on_spot$spot)
    cell_id = as.character(cells_on_spot$cell_id)
    
    PE_on_cell = PE_on_spot[ match(spot_id,rownames(PE_on_spot)) ,  ]
    PME_on_cell = t( apply(as.matrix(PM_on_cell*PE_on_cell),1,function(x)x/sum(x,na.rm=T)) )
    
    PME_uni_cell = apply(as.matrix(PME_on_cell), 2, function(x) tapply(x, cell_id, max) )
    PME_uni_cell = PME_uni_cell[ match( cell_id, rownames(PME_uni_cell) ) , ]
    
    a = t(PME_uni_cell) %*% as.matrix(cells_on_spot[,features])
    s = colSums(PME_uni_cell, na.rm=T)
    mu = a/s
    
    b = do.call( rbind, lapply(1:ncol(PM_on_cell), function(t) {
        p = as.matrix(PME_uni_cell[,t])
        if(length(features)>1) d = t(apply(cells_on_spot[,features],1,function(x) x-mu[t,] ))
        if(length(features)==1) d = matrix( cells_on_spot[,features]-mu[t,], ncol=1)
        t(p)%*%d^2
    }))
    rownames(b) = colnames(PM_on_cell)
    #sigma = apply(b,2,function(x) x/s)
    sigma = sqrt( b/s )
    
    colnames(mu) = colnames(sigma) = features
    list(mu=mu, sigma=sigma)
    
}
