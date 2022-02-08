#' STIE
#'
#' @param ST_expr 
#' @param Signature 
#' @param cells_on_spot 
#' @param features 
#' @param lambda 
#' @param steps 
#' @param morphology_steps 
#' @param known_signature 
#' @param known_cell_types 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
STIE <- function(ST_expr, Signature, cells_on_spot, features, 
                  lambda=0, steps=30, morphology_steps=ceiling(steps/3),
                  known_signature=TRUE, known_cell_types=FALSE)
{
    # S is the matrix of cell type signature, row: genes, col: cell types
    # B is the gene expression of one spot, row: genes, col: cell types
    # beta is the current estimation of regression coefficients
    # Pm is the probability of morphology for cell type
    # lambda is the langurange multiplier
    
    cl <- makeCluster( detectCores() )
    registerDoParallel(cl)
    
    #cells_on_spot <- get_cells_on_spot( cell_coordinates=cell_coordinates, spot_coordinates, 2*spot_radius)
    spot_id = as.character(cells_on_spot$spot)
    cell_id = as.character(cells_on_spot$cell_id)
    
    ST_expr2 = t( ST_expr[ rownames(ST_expr)%in%spot_id, match(rownames(Signature),colnames(ST_expr)) ] )
    
    ##########################################################################################
    ########### set initial values
    ##########################################################################################
    if( known_signature & !known_cell_types ) 
    {
        Signature_ = Signature
        
        ST_expr2 = t( ST_expr[ rownames(ST_expr)%in%spot_id, match(rownames(Signature),colnames(ST_expr)) ] )
        
        mu = sigma = matrix(nrow=ncol(Signature),ncol=length(features),data=0)
        rownames(mu) = rownames(sigma) = colnames(Signature)
        colnames(mu) = colnames(sigma) = features
        
        PE_on_spot = t( apply( ST_expr2, 2, function(x) solveOLS( Signature, as.matrix(x) ) ) )
        
        for(t in 1:nrow(mu))
        {
            ford = order(PE_on_spot[,t],decreasing=T)
            mu[t,] = apply( subset(cells_on_spot,spot%in%rownames(PE_on_spot)[ford[1:40]])[,features], 2, mean ) 
            sigma[t,] = apply( subset(cells_on_spot,spot%in%rownames(PE_on_spot)[ford[1:40]])[,features], 2, sd ) 
        }
        
        PM_on_cell = calculate_morphology_probability(cells_on_spot, features, mu, sigma )
        
    }
    
    if( !known_signature & !known_cell_types ) 
    {
        Signature_ = Signature
        
        ST_expr2 = t( ST_expr[ rownames(ST_expr)%in%spot_id, match(rownames(Signature),colnames(ST_expr)) ] )
        
        mu = sigma = matrix(nrow=ncol(Signature),ncol=length(features),data=0)
        rownames(mu) = rownames(sigma) = colnames(Signature)
        colnames(mu) = colnames(sigma) = features
        
        PE_on_spot = t( apply( ST_expr2, 2, function(x) solveOLS( Signature, as.matrix(x) ) ) )
        
        for(t in 1:nrow(mu))
        {
            ford = order(PE_on_spot[,t],decreasing=T)
            mu[t,] = apply( subset(cells_on_spot,spot%in%rownames(PE_on_spot)[ford[1:40]])[,features], 2, mean ) 
            sigma[t,] = apply( subset(cells_on_spot,spot%in%rownames(PE_on_spot)[ford[1:40]])[,features], 2, sd ) 
        }
        
        PM_on_cell = calculate_morphology_probability(cells_on_spot, features, mu, sigma )
        
    }
    
    if( !known_signature & known_cell_types ) 
    {
        mu = do.call(cbind, lapply(features, function(f) tapply( cell_coordinates[,f], 
                                                  as.character(cell_coordinates[,"celltypes"]), mean ) ))
        sigma = do.call(cbind, lapply(features, function(f) tapply( cell_coordinates[,f], 
                                                                 as.character(cell_coordinates[,"celltypes"]), sd ) ))
        colnames(mu) = colnames(sigma) = features
        
        PM_on_cell = calculate_morphology_probability(cells_on_spot, features, mu, sigma )
        coefs = apply(PM_on_cell, 2, function(x) tapply(x,spot_id,sum) )
        Signature = t( apply( ST_expr2, 1, function(x) solveOLS( coefs, as.matrix(x), scaled=F ) ) )
        
        PE_on_spot = t( apply(coefs, 1, function(x) x/sum(x)) )
    }
    
    ##########################################################################################
    
    step = 0
    t1 = Sys.time()
    mu_ = mu
    sigma_ = sigma
    Signature_ = Signature
    PE_on_spot_ = PE_on_spot
    
    while(step<steps) {
        
        step = step+1
        t0 = t1
        cat( "step", step, "##########\n" )
        
        ########### 
        cat( "   recalculate PE and PM ... \n")
        mu = mu_
        sigma = sigma_
        Signature = Signature_
        
        PE_on_spot = PE_on_spot_
        PM_on_cell = calculate_morphology_probability(cells_on_spot, features, mu, sigma )
        
        ########### 
        PE_on_cell = PE_on_spot[ match(spot_id,rownames(PE_on_spot)) ,  ]
        
        # v1.4
        #PME_on_cell = t( apply(PM_on_cell*PE_on_cell,1,function(x)x/sum(x)) )
        #PM_on_spot = apply(PME_on_cell, 2, function(x) tapply(x,spot_id,sum) )
        
        # v1.5
        PM_on_cell = t( apply(PM_on_cell, 1, function(x)x/sum(x)) )
        PM_on_spot = apply(PM_on_cell, 2, function(x) tapply(x,spot_id,sum) )
        
        PM_on_spot = PM_on_spot[ match(rownames(PE_on_spot),rownames(PM_on_spot)), ]
        PM_on_spot = t( apply( PM_on_spot, 1, function(x) x/sum(x) ) )
        all( rownames(PE_on_spot)==rownames(PM_on_spot) )
        
        ###########  
        if(step>10)
        {
            cat( "   updating morphology parameter ... \n")
            morphology_parameter <- update_morphology_parameter(PE_on_spot, PM_on_cell, cells_on_spot, features)
            mu_ = morphology_parameter$mu
            sigma_ = morphology_parameter$sigma
        }
        ########### 
        cat( "   updating expression regression parameter ... \n" ) 
        coefs = do.call(rbind, foreach( i=1:nrow(PE_on_spot), .packages=c("quadprog","STIE")) %dopar% {
            
            Expr_on_spot_i = as.matrix( ST_expr2[,i] )
            PE_on_spot_i = PE_on_spot[i,]
            PM_on_spot_i = PM_on_spot[i,]
            update_expression_regression_parameter(Signature, 
                                                   Expr_on_spot_i, PE_on_spot_i, 
                                                   PM_on_spot_i, lambda = lambda, 
                                                   scaled=F, kfold=10 )
        })
        
        PE_on_spot_ = t( apply(coefs, 1, function(x) x/sum(x)) )
        rownames(coefs) = rownames(PE_on_spot)
        rownames(PE_on_spot_) = rownames(PE_on_spot)
        
        if( !known_signature )
        {
            Signature_ = t( apply( ST_expr2, 1, function(x) solveOLS( coefs, as.matrix(x), scaled=F ) ) )
            #plot(Signature_[,1],Signature[,1])
        }
        
        Sig_tah = sum( abs(Signature_ - Signature) )/nrow(Signature_)
        
        PE_tah = sum( abs(PE_on_spot_ - PE_on_spot) )/nrow(PE_on_spot)
        PM_tah = sum( abs(mu_ - mu) )/nrow(mu)
        
        t1 = Sys.time()
        print(t1-t0)
        cat("--> PE diff:", PE_tah, ", and PM diff:", PM_tah, ", and Sig diff:", Sig_tah, "\n\n")
        
        print(mu)
    }
    

    stopCluster(cl)
    
    PME_on_cell = t( apply(PM_on_cell*PE_on_cell,1,function(x)x/sum(x)) )
    PME_uni_cell = apply(PME_on_cell, 2, function(x) tapply(x, cell_id, max) )
    PME_uni_cell = PME_uni_cell[ match( cell_id, rownames(PME_uni_cell) ) , ]
    cell_types = colnames(PME_uni_cell)[apply(PME_uni_cell, 1, which.max)]
    names(cell_types) = rownames(PME_uni_cell)
    
    cat( "   calculate adjacent similarity ... \n")
    #adjacent_similarity = get_adjacent_similarity(cells_on_spot, features, cell_types) 
    
    result <- list(lambda=lambda,
         mu=mu, 
         sigma=sigma, 
         PE_on_spot = PE_on_spot,
         PM_on_cell = PM_on_cell,
         PME_uni_cell = PME_uni_cell, 
         cell_types=cell_types, 
         Signature=Signature,
         cells_on_spot=cells_on_spot)
    
}
