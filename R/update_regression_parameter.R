#' update_expression_regression_parameter
#'
#' update_expression_regression_parameter re-estimates the parameters for the gene expression model in STIE
#'
#' @param Signature a matrix of gene expression in gene X cell type, representing cell type transcriptomic signatures 
#' @param Expr_on_spot_i a vector of numeric values representing the gene expression on the i-th spot
#' @param PE_on_spot_i a matrix of numeric values representing the probability of gene expression on the i-th spot
#' @param PM_on_spot_i a matrix of numeric values representing the probability of nuclear morphology on the i-th spot
#' @param lambda a numeric value representing shrinkage penalty of nuclear morphology
#' @param scaled a boolean value representing whether to scale the probability to be sum of 1
#' @param kfold an integer value for K-fold cross validation 
#'
#' @return a vector of numeric values, representing cell type proportion of the i-th spot
#' @export
#'
#' @author Shijia Zhu, \email{shijia.zhu@@utsouthwestern.edu}
#'
#' @references
#'
#' @seealso \code{\link{solveNNLS}}; 
#' 
#' 
update_expression_regression_parameter <- function( Signature, Expr_on_spot_i, 
                                                    PE_on_spot_i, PM_on_spot_i, 
                                                    lambda = c(0,10^c(0:7)), scaled=T, 
                                                    kfold=10 ) {
    
    
    solveOLS2 <- function(Signature, Expr_on_spot_i, 
                          PE_on_spot_i, PM_on_spot_i, 
                          lambda, scaled=T) {
        
        # Signature is the matrix of cell type signature
        # Expr_on_spot_i is the gene expression of i-th spot
        # PE_on_spot_i is the current estimation of regression coefficients
        # PM_on_spot_i is the probability of morphology for cell type
        # lambda is the langurange multiplier
        
        if(0) {
            s2 = var( as.numeric(Expr_on_spot_i) - Signature%*%as.numeric(PE_on_spot_i) )
            s2 = 1
            lambda_ = as.numeric( lambda*s2/sum(PE_on_spot_i)^2 )
            cat("lambda_:", lambda_, s2, "\n")
        }
        
        lambda_ = as.numeric( lambda/sum(PE_on_spot_i)^2 )
        
        t <- ncol(Signature)
        I <- diag(t)
        
        D <- t(Signature)%*%Signature + lambda_*I
        d <- t(Signature)%*%Expr_on_spot_i + lambda_*sum(PE_on_spot_i)*I%*%PM_on_spot_i
        A <- I
        bzero <- c(rep(0,t))
        solution <- solve.QP(D,d,A)$solution
        names(solution) <- colnames(Signature)
        
        solution[solution<0] = 1e-100
        
        if(scaled) solution <- solution/sum(solution)
        solution
    }
    
    
    if( length(lambda)==1 ) {
        
        solution = solveOLS2(Signature, Expr_on_spot_i, PE_on_spot_i, PM_on_spot_i, lambda, scaled=scaled)
        
    } else {
        
        set.seed(1234)
        
        n = nrow(Expr_on_spot_i)
        m = ceiling( n/kfold )
        index = sample(1:n,n)
        
        sst = rep(0,length(lambda))
        sse = rep(0,length(lambda))
        
        for(j in 1:length(lambda))
        {
            #cat(j,"\n")
            for(k in 1:kfold) {
                
                s = (k-1)*m+1
                e = min( k*m, n)
                index_test = index[s:e]
                index_train = setdiff(index, index_test)
                
                solution = solveOLS2(Signature[index_train,], as.matrix(Expr_on_spot_i[index_train]), 
                                     PE_on_spot_i, PM_on_spot_i, 
                                     lambda[j], scaled=F)
                predicted = Signature[index_test,]%*%solution
                
                sst[j] = sst[j] + sum( ( Expr_on_spot_i[index_test] - mean(Expr_on_spot_i[index_test]) )^2 )
                sse[j] = sse[j] + sum( ( predicted - Expr_on_spot_i[index_test] )^2 )
                
            }
            
        }
        
        # rsq <- 1 - sse/sst
        
        j = which.min(sse)
        solution = solveOLS2(Signature, Expr_on_spot_i, PE_on_spot_i, PM_on_spot_i, lambda[j], scaled=scaled)
    }
    
    solution
    
}


