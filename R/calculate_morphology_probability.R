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

