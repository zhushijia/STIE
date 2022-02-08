#' update_expression_regression_parameter
#'
#' @param Signature 
#' @param Expr_on_spot_i 
#' @param PE_on_spot_i 
#' @param PM_on_spot_i 
#' @param lambda 
#' @param scaled 
#' @param kfold 
#'
#' @return
#' @export
#'
#' @examples
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
        
        t <- ncol(Signature)
        I <- diag(t)
        
        D <- t(Signature)%*%Signature + lambda*I
        d <- t(Signature)%*%Expr_on_spot_i + lambda*sum(PE_on_spot_i)*I%*%PM_on_spot_i
        A <- I
        bzero <- c(rep(0,t))
        solution <- solve.QP(D,d,A)$solution
        names(solution) <- colnames(Signature)
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


