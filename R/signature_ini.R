#' signature_ini assign initial values of transcriptomic signature for STIE clustering
#'
#' @param k integer value indicating the number of clusters
#' @param ST_expr a matrix of gene expression representing the spot level gene expression, with row being the spot and column being the gene
#' @param cells_on_spot a data frame representing the cells on spots along with the cellular morphological features. It will be used when on=='morph'
#' @param on a character value taking from c("expr",'morph'), indicating whether to use gene expression initial value at spot level or morphological value at the single-cell level
#' @param features a vector of character values indicating what mophorlogical features to use, when on=='morph', 
#' @param pc_n an integer value indicating the number of PCs used for K-means on gene expression, when on=='expr' (by default pc_n=10) 
#' @param rho_cutoff a numeric value indicating the cutoff for the similarity between transcriptomic signatures. Once greater than the cutoff, one of the signature is removed, due to redundancy.  
#'
#' @return a matrix of gene expression representing the cell type transcriptomic signature, with row being the gene and column being the cell type
#' @export
#'
#' @examples
#' Signature_ = signature_ini(k=5, ST_expr=ST_expr, on="expr")
#' Signature_ = signature_ini(k=5, cells_on_spot=cells_on_spot, on='morph', features=c('size','shape'))
#' 
#' 
signature_ini = function(k, ST_expr=NULL, cells_on_spot=NULL, on=c("expr",'morph'),
                         features=c('size','shape'), pc_n = 10, rho_cutoff=0.99 ) 
{
    on = match.arg(on)
    
    if(on=="expr") {
        pc = prcomp(ST_expr)$x[,1:pc_n]
        set.seed(1234)
        cluster = kmeans(pc,k)$cluster
        Signature_ = t(apply(ST_expr, 2, function(x) tapply(x,cluster,mean) ))
    }
    
    if(on=="morph") {
        set.seed(1234)
        km = kmeans(cells_on_spot[,features], k)$cluster
        coefs = table( cells_on_spot$spot, km )
        coefs = coefs[rownames(ST_expr), ]
        Signature_ = t( apply( ST_expr, 2, function(x) solveNNLS( coefs, as.matrix(x), scaled=F ) ) )
    }
    
    cc = cor(Signature_)
    diag(cc) = 0
    while( any(cc>rho_cutoff) ) {
        ind = which(cc>rho_cutoff, arr.ind=T)
        j = ind[1,1]
        cc = cc[-j, -j]
    }
    
    Signature_[,colnames(cc)]
    
}



