#' solveOLS
#'
#' @param S 
#' @param B 
#' @param scaled 
#'
#' @return
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

