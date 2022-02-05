#' update_morphology_parameter
#'
#' @param PE_on_spot 
#' @param PM_on_cell 
#' @param cells_on_spot 
#' @param features 
#'
#' @return
#' @export
#'
#' @examples
update_morphology_parameter <- function(PE_on_spot, PM_on_cell, cells_on_spot, features)
{
    spot_id = as.character(cells_on_spot$spot)
    cell_id = as.character(cells_on_spot$cell_id)
    
    PE_on_cell = PE_on_spot[ match(spot_id,rownames(PE_on_spot)) ,  ]
    PME_on_cell = t( apply(PM_on_cell*PE_on_cell,1,function(x)x/sum(x)) )
    
    PME_uni_cell = apply(PME_on_cell, 2, function(x) tapply(x, cell_id, max) )
    PME_uni_cell = PME_uni_cell[ match( cell_id, rownames(PME_uni_cell) ) , ]
    
    a = t(PME_uni_cell) %*% as.matrix(cells_on_spot[,features])
    s = colSums(PME_uni_cell)
    mu = a/s
    
    b = do.call( rbind, lapply(1:ncol(Signature), function(t) {
        p = as.matrix(PME_uni_cell[,t])
        d = t(apply(cells_on_spot[,features],1,function(x) x-mu[t,] ))
        t(p)%*%d^2
    }))
    rownames(b) = colnames(Signature)
    #sigma = apply(b,2,function(x) x/s)
    sigma = sqrt( b/s )
    
    list(mu=mu, sigma=sigma)
    
}
