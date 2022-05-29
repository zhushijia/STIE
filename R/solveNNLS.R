#' solveNNLS
#'
#' solveNNLS solves the linear regression model with non-negative regression coefficients
#'
#' @param S a matrix of cell type transcriptomic signature with row representing the gene and column representing the cell type 
#' @param B a vector of bulk gene expression with row representing the gene and column representing the sample
#' @param scaled a boolean value indicating whether to scale the sum of cell type proportions to be 1
#'
#' @return a vector of non-negative regression coefficients on cell type transcriptomic signatures
#' @export
#'
#' @examples
#' 
#' 
solveNNLS <- function(S, B, scaled=T) {
    D<-t(S)%*%S
    d<-t(S)%*%B
    A<-cbind(diag(dim(S)[2]))
    bzero<-c(rep(0,dim(S)[2]))
    solution<-solve.QP(D,d,A,bzero)$solution
    names(solution)<-colnames(S)
    solution = abs(solution)
    if(scaled) solution = solution/sum(solution)
    return(solution)
}

