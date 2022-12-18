#' STIE_search
#'
#' STIE_search: hyperparameter and feature selection using STIE function
#'
#' @param ST_expr a matrix of gene expression representing the spot level gene expression, with row being the spot and column being the gene
#' @param Signature a matrix of gene expression representing the cell type transcriptomic signature, with row being the gene and column being the cell type
#' @param cells_on_spot a data frame representing the cells on spots along with the cellular morphological features, where "cell_id" represents the unique cell id and "spot" represents the uniqe spot id
#' @param features a vector of character values, representing the morphological features used in the STIE model
#' @param steps an integer value representing the number of iterative steps of EM algorithm 
#' @param morphology_steps an integer value representing the number of iterative steps of updating mophological model before updating gene expression model
#' @param known_signature a boolean value representing whether the transcriptomic signature is given. If TRUE, STIE performs single cell deconvolution
#' @param known_cell_types a boolean value representing whether the cell typing is given. If TRUE, STIE will perform cell type signature learning
#' @param min_cells  a boolean value representing the minimum number of cells to keep for each cell type. If the cell count is smaller than min_cells, the cell type will be eliminated
#' @param lambdas a vector of numeric value representing the shrinkage penalty of nuclear morphology
#' @param equal_prior  a boolean value representing the wheter or not assume equal prior, i.e., P(q|theta)
#' @param plot a boolean value representing whether to plot the result
#' 
#'
#' @return A list containing the follow components:
#' \itemize{
#'  \item {lambda} a numeric value representing the shrinkage penalty of nuclear morphology, which is the same data frame cells_on_spot with the input
#'  \item {mu} a vector of numeric values representing the means of morphological features
#'  \item {sigma} a vector of numeric values representing the standard deviations of morphological features
#'  \item {PE_on_spot} a matrix of cell type probabilities for each spot estimated from the spot gene expression, where the row is the spot and column is the cell type
#'  \item {PM_on_cell} a matrix of cell type probabilities for each cell estimated from the cellular morphological features, where the row is the cell and column is the cell type
#'  \item {PME_uni_cell} a matrix of cell type probabilities for each cell estimated from both spot gene expression and cellular morphological features, where the row is the cell and column is the cell type
#'  \item {cell_types} a vector of character values representing the cell type for each cell, where the names(cell_types) are the unique cell ids
#'  \item {uni_cell_types} a vector of character values representing the cell type of unqiue cells with redundant cells eliminated 
#'  \item {Signature} a matrix of gene expression in gene X cell type. During deconvolution, it is the same with the input. During clustering, it is the re-estimated signature from the ST data. 
#'  \item {cells_on_spot} a data frame representing the cells on spots along with the cellular morphological features, which is the same data frame cells_on_spot with the input
#' }
#' 
#' @author Shijia Zhu, \email{shijia.zhu@@utsouthwestern.edu}
#' @export
#'
#' @references
#'
#' @seealso \code{\link{get_cells_on_spot}}; \code{\link{split_image}}; \code{\link{run_imageJ_plugin}}; \code{\link{merge_feature}};
#' 
#' 
STIE_search0 <- function(ST_expr, Signature, cells_on_spot, 
                        steps=30, morphology_steps=ceiling(steps/3),
                        known_signature=TRUE, known_cell_types=FALSE, min_cells=2, equal_prior=TRUE,
                        lambdas=c(0,1e3,1e6), features_of_interest = c("size", "shape"),
                        criterion = c("L2sum","rmse", "logLik") , 
                        plot=TRUE) {
    
    criterion = match.arg(criterion)
    
    features_list <- list( size = c("Area", 'Major', 'Minor', 'Width', 'Height', 'Feret','Perim.'),
                           shape = c("Round", 'Circ.'), 
                           angle = c('FeretAngle','Angle'),
                           solidity = c("Solidity") )
    
    if(any( names(features_list) %in% features_of_interest )) {
        PCs <- do.call(cbind, lapply(features_list, function(f) {
            X = cells_on_spot[, f]
            X2 = scale(X)
            prcomp(X2)$x[, 1]
        }))
        
        cells_on_spot2 <- cbind(cells_on_spot, PCs)
        
    } else {
        
        cells_on_spot2 <- cells_on_spot
    }
    
    features_of_interest = intersect(features_of_interest, colnames(cells_on_spot2))
    stopifnot( length(features_of_interest)>0  )
    
    paths <- list()
    
    for( i in 1:length(lambdas) ) {
        
        lai = lambdas[i]
        
        result1 = lapply( features_of_interest, function(cf) 
            STIE(ST_expr, Signature, cells_on_spot2, features=cf, 
                 lambda=lai, steps=steps, 
                 known_signature=known_signature, known_cell_types=known_cell_types,
                 min_cells=min_cells, equal_prior=equal_prior) )
        score1 = lapply(result1, function(result) get_summary(result,ST_expr) )
        names(score1) = names(result1) = features_of_interest
        
        if(criterion=='rmse')   ord1 = order( sapply(score1,function(x)x$rmse) )
        if(criterion=='L2sum')     ord1 = order( sapply(score1,function(x)x$L2sum) )
        if(criterion=='logLik') ord1 = order( sapply(score1,function(x)x$logLik), decreasing=T )
        
        result1 = result1[ord1]
        score1 = score1[ord1]
        
        ################################################
        results = list()
        scores = list()
        results[[1]] = result1[[1]]
        scores[[1]] = score1[[1]]
        names(results)[1] = names(scores)[1] = features_of_interest[ord1[1]]
        
        if(length(features_of_interest)>1) {
            for(j in 2:length(features_of_interest)) {
                results[[j]] = STIE(ST_expr, Signature, cells_on_spot2, features=features_of_interest[ord1[1:j]], 
                                    lambda=lai, steps=steps, 
                                    known_signature=known_signature, known_cell_types=known_cell_types,
                                    min_cells=min_cells, equal_prior=equal_prior )
                scores[[j]] = get_summary(results[[j]], ST_expr)
                names(results)[j] = names(scores)[j] = paste(features_of_interest[ord1[1:j]],collapse="+")
            }
        }
       
        if(0) {
            paths[[i]] = list( ordered_features=names(result1), 
                               results_ord1 = result1, scores_ord1 = score1,
                               results_ord2 = results, scores_ord2 = scores,
                               val_ord1 = sapply(score1,function(x)x[[criterion]]),
                               val_ord2 = sapply(scores,function(x)x[[criterion]]),
                               logLik_ord1 = sapply(score1,function(x)x$logLik),
                               logLik_ord2 = sapply(scores,function(x)x$logLik),
                               logLikMorph_ord1 = sapply(score1,function(x) mean(x$logLik_Morp) ),
                               logLikMorph_ord2 = sapply(scores,function(x) mean(x$logLik_Morp) ),
                               criterion=criterion )
        }
        
        paths[[i]] = list( ordered_features=names(result1), 
                           results = results, 
                           scores = scores,
                           values = sapply(scores,function(x)x[[criterion]]),
                           logLik = sapply(scores,function(x)x$logLik),
                           logLikMorph = sapply(scores,function(x) mean(x$logLik_Morp) ),
                           criterion=criterion )
        
        names(paths)[i] = lai
        
    }
    
    if(plot) {
        par(mfrow=c(2,2))
        par(mar=c(2,2,2,2))
        plot_searchPath(paths)
        #plot_pathScore( lapply(paths,function(x)x$logLik), name='logLik')
        plot_pathScore( lapply(paths,function(x)x$logLikMorph), name='logLik_Morph')
        boxplot( lapply(paths, function(x) sapply(x$scores, function(y) y$L2sum) ), main='L2sum')
        boxplot( lapply(paths, function(x) x$logLikMorph), main='logLik_Morph')
    }
    
    paths
    
}


plot_searchPath <- function(paths)
{
    #par(mfrow=c(1,1))
    par(mar=c(2,2,2,2))
    best_val = lapply(paths, function(x) x$values)
    n = length(paths)
    f = length(best_val[[1]])
    ylim = range(do.call(c,best_val))
    ylim = c( ylim[1], ylim[1] + 1.5*(ylim[2]-ylim[1]) )
    plot(NA, xlim=c(0.5,f+1), ylim=ylim, 
         xlab="# Features", ylab='RMSE', main=paste0( paths[[1]]$criterion, ': Features & Lambda'), 
         axes = FALSE) #xaxt="n", yaxt="n" )
    axis(1, at = c(1:f))
    axis(2)
    for( i in 1:n ) lines(best_val[[i]],col=i)
    for( i in 1:n ) points(best_val[[i]])
    for( i in 1:n ) points(which.min(best_val[[i]]), min(best_val[[i]]), col=i, pch=8, cex=1.5 )
    for( i in 1:n ) text(1:f,best_val[[i]],names(best_val[[i]]),pos=3)
    legend( 'topright', legend=paste0('lambda=',names(best_val)), col=1:n, lty=1, bty = "n")
}

plot_pathScore <- function(score,name=NULL)
{
    par(mar=c(2,2,2,2))
    n = length(score)
    f = length(score[[1]])
    ylim = range(do.call(c,score))
    ylim = c( ylim[1], ylim[1] + 1.5*(ylim[2]-ylim[1]) )
    plot(NA, xlim=c(0.5,f+1), ylim=ylim, 
         xlab="# Features", ylab=name, main=paste0(name,': Features & Lambda'), 
         axes = FALSE) #xaxt="n", yaxt="n" )
    axis(1, at = c(1:f))
    axis(2)
    for( i in 1:n ) lines(score[[i]],col=i)
    for( i in 1:n ) points(score[[i]])
    #for( i in 1:n ) points(which.min(score[[i]]), min(score[[i]]), col=i, pch=8, cex=1.5 )
    for( i in 1:n ) text(1:f,score[[i]],names(score[[i]]),pos=3)
    legend( 'topright', legend=paste0('lambda=',names(score)), col=1:n, lty=1, bty = "n")
}



