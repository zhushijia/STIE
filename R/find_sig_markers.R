#' find_signature_markers
#'
#' @param STIE_result 
#' @param database 
#' @param db_category 
#' @param DEG_pthres 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
#' 
find_sig_markers <- function(STIE_result, ST_expr, transform="inverse.transform", DEG_pthres=1e-5)
{
    # https://github.com/sqjin/CellChat
    # http://www.cellchat.org/
    
    Signature = STIE_result$Signature
    cells_on_spot = STIE_result$cells_on_spot
    PE_on_spot = STIE_result$PE_on_spot
    PM_on_cell = STIE_result$PM_on_cell
    cell_types = STIE_result$cell_types
    lambda = STIE_result$lambda
    
    spot_id = as.character(cells_on_spot$spot)
    cell_id = as.character(cells_on_spot$cell_id)
    PM_on_spot = apply(PM_on_cell, 2, function(x) tapply(x,spot_id,sum) )
    PM_on_spot = PM_on_spot[ match(rownames(PE_on_spot),rownames(PM_on_spot)), ]
    PM_on_spot = t( apply( PM_on_spot, 1, function(x) x/sum(x) ) )
    all( rownames(PE_on_spot)==rownames(PM_on_spot) )
    
    ST_expr2 = t( ST_expr[ rownames(ST_expr)%in%spot_id, match(rownames(Signature),colnames(ST_expr)) ] )
    ST_expr3 = t( ST_expr[ rownames(ST_expr)%in%spot_id, ] )
    
    coefs = do.call(rbind, lapply( 1:nrow(PE_on_spot), function(i) {
        
        #cat(i,"\n")
        Expr_on_spot_i = as.matrix( ST_expr2[,i] )
        PE_on_spot_i = PE_on_spot[i,]
        PM_on_spot_i = PM_on_spot[i,]
        update_expression_regression_parameter(Signature, 
                                               Expr_on_spot_i, PE_on_spot_i, 
                                               PM_on_spot_i, lambda = lambda, 
                                               scaled=F, kfold=10 )
    }))
    
    rownames(coefs) = rownames(PE_on_spot)
    colnames(coefs) = colnames(Signature)
    
    Signature3 = t( apply( ST_expr3, 1, function(x) solveNNLS( coefs, as.matrix(x), scaled=F ) ) )
    
    INT = function(x, k=3/8)
    {
        x = as.numeric(x)
        a = rank(x,na.last="keep") - k
        b = sum(!is.na(x)) - 2*k + 1
        qnorm( a/b )
        #qnorm( ( rank(x,na.last="keep")-0.5 ) / sum(!is.na(x)) )
    }
    
    ttest.2sample <- function( m1, sd1, n1, m2, sd2, n2 ) {
        m = m1-m2
        se = sqrt( sd1^2/n1 + sd2^2/n2 )
        t = m/se
        df = ( sd1^2/n1 + sd2^2/n2 )^2 / ( (sd1/n1)^2/(n1-1) + (sd2/n2)^2/(n2-1) )
        p = pt( t, df=df, lower.tail=F )
        # cat(m,se,t,df,p,"\n")
        data.frame( m=m, m1=m1, m2=m2, t=t, df=df, p=p )
    }
    
    features.info = do.call(rbind, lapply( 1:nrow(ST_expr3), function(i) {
        
        cat(i,"\n")
        
        y = ST_expr3[i,]
        if(transform=="log") y = log(y+1)
        if(transform=="inverse.transform") y = INT( y )

        z = summary(lm( y ~ coefs-1 ))$coef
        n = colSums(coefs)
        
        if(0)
        {
            m0 = mean(y)
            sd0 = sd(y)
            n2=length(y)
        } 
        
        info = do.call(rbind, lapply( 1:ncol(Signature), function(j) {
            
            m0 = sum(z[-j,1]*n[-j])/sum(n[-j])
            sd0 = sqrt( sum( z[-j,2]^2*n[-j] ) )
            n2 = sum(n[-j])
            
            with( ttest.2sample( m1=z[j,1], sd1=sqrt( z[j,2]^2*n[j] ), n1=n[j], 
                                m2=m0, sd2=sd0, n2=n2 ), 
                 data.frame( clusters = colnames(Signature)[j],
                             features = rownames(ST_expr3)[i],
                             diff = m,
                             pct.1 = m1,
                             pct.2 = m2,
                             pvalues = p))
        }))
        
        subset(info, pvalues<DEG_pthres)
        
    } ))
    
    features = unique(as.character(features.info$features))
    var.features = list(features.info=features.info,
                        features = features)
    
    features.info
    
}

