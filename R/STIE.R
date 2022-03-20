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
#' @param min_cells 
#' 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
STIE <- function(ST_expr, Signature, cells_on_spot, features, 
                  lambda=0, steps=30, morphology_steps=ceiling(steps/3),
                  known_signature=TRUE, known_cell_types=FALSE, min_cells=2)
{
    
    cl <- makeCluster( detectCores() )
    registerDoParallel(cl)
    
    #cells_on_spot <- get_cells_on_spot( cell_coordinates=cell_coordinates, spot_coordinates, 2*spot_radius)
    spot_id = as.character(cells_on_spot$spot)
    cell_id = as.character(cells_on_spot$cell_id)
    
    if( is.null(Signature) )
    {
        ST_expr2 = t( ST_expr[ rownames(ST_expr)%in%spot_id,  ] )
    } else {
        Signature = Signature[rownames(Signature)%in%colnames(ST_expr),]
        ST_expr2 = t( ST_expr[ rownames(ST_expr)%in%spot_id, match(rownames(Signature),colnames(ST_expr)) ] )
    }
    
    ##########################################################################################
    ########### set initial values
    ##########################################################################################
    if( known_signature & !known_cell_types )  # deconvolution
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
    
    if( !known_signature & !known_cell_types ) # clustering
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
        annoted_cells = which(!is.na(cells_on_spot$cell_types))
        mu = do.call(cbind, lapply(features, function(f) tapply( cells_on_spot[,f], 
                                                  as.character(cells_on_spot[,"cell_types"]), mean ) ))
        sigma = do.call(cbind, lapply(features, function(f) tapply( cells_on_spot[,f], 
                                                                 as.character(cells_on_spot[,"cell_types"]), sd ) ))
        colnames(mu) = colnames(sigma) = features
        
        coefs = table( as.character(cells_on_spot$spot), as.character(cells_on_spot$cell_types) )
        coefs = coefs[match( colnames(ST_expr2), rownames(coefs) ), ]
        Signature = t( apply( ST_expr2, 1, function(x) solveOLS( coefs, as.matrix(x), scaled=F ) ) )
        
        PM_on_cell = calculate_morphology_probability(cells_on_spot, features, mu, sigma )
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
        PE_on_cell = PE_on_spot[ match(spot_id,rownames(PE_on_spot)) ,  ]
        
        ####################################################################################
        ###########  v1.8 to filter out cell types of small number
        ####################################################################################
        #PME_on_cell = t( apply(PM_on_cell*PE_on_cell,1,function(x)x/sum(x)) )
        PME_on_cell = PM_on_cell*PE_on_cell
        PME_on_cell = do.call(rbind, lapply( 1:nrow(PME_on_cell), function(i) {
            s = sum(PME_on_cell[i,])
            if( s>0 ) {
                PME_on_cell[i,]/s
            } else {
                PE_on_cell[i,]
            } } ) )
        PME_uni_cell = apply(PME_on_cell, 2, function(x) tapply(x, cell_id, max) )
        PME_uni_cell = PME_uni_cell[ match( cell_id, rownames(PME_uni_cell) ) , ]
        cell_types = colnames(PME_uni_cell)[apply(PME_uni_cell, 1, which.max)]
        names(cell_types) = rownames(PME_uni_cell)
        
        uni_cell_types = tapply( cell_types, names(cell_types), function(x) x[1] )
        count_cell_types = table(uni_cell_types)
        filtered = which( colnames(Signature) %in% names(count_cell_types)[count_cell_types > min_cells] )
        
        if( length(filtered)<ncol(Signature) )
        {
            mu = mu_[filtered, ]
            sigma = sigma_[filtered, ]
            Signature = Signature_[ ,filtered]
            PE_on_spot = t( apply(PE_on_spot_[ ,filtered], 1, function(x) x/sum(x)) )
            PM_on_cell = calculate_morphology_probability(cells_on_spot, features, mu, sigma )
            PE_on_cell = PE_on_spot[ match(spot_id,rownames(PE_on_spot)) ,  ]
            
            mu_ = mu
            sigma_ = sigma
            Signature_ = Signature
            PE_on_spot_ = PE_on_spot
            
        } 
        
        ####################################################################################
        ###########  above v1.8 to filter out cell types of small number
        ####################################################################################
        
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
        if(step>morphology_steps)
        {
            cat( "   updating morphology parameter ... \n")
            morphology_parameter <- update_morphology_parameter(PE_on_spot, PM_on_cell, cells_on_spot, features)
            mu_ = morphology_parameter$mu
            sigma_ = morphology_parameter$sigma
        }
        ########### 
        cat( "   updating expression regression parameter ... \n" ) 
        #ccoefs = do.call(rbind, foreach( i=1:nrow(PE_on_spot), .packages=c("quadprog","STIE")) %dopar% {
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
        
        PE_on_spot_ = t( apply(coefs, 1, function(x) x/sum(x)) )
        rownames(coefs) = rownames(PE_on_spot)
        rownames(PE_on_spot_) = rownames(PE_on_spot)
        
        if( !known_signature )
        {
            # Signature_ = t( apply( ST_expr2, 1, function(x) solveOLS( coefs, as.matrix(x), scaled=F ) ) )
            Signature_ = do.call( rbind, lapply( 1:nrow(ST_expr2), function(i) {
                tryCatch( solveOLS( coefs, as.matrix(ST_expr2[i,]), scaled=F ), 
                          error = function(e) {
                              rep(0,ncol(coefs))
                          })
            }))
            rownames(Signature_) = rownames(ST_expr2)
            #plot(Signature_[,1],Signature[,1])
        }
        
        Sig_tah = sum( abs(Signature_ - Signature) )/nrow(Signature_)
        
        PE_tah = sum( abs(PE_on_spot_ - PE_on_spot) )/nrow(PE_on_spot)
        PM_tah = sum( abs(mu_ - mu) )/nrow(mu)
        
        t1 = Sys.time()
        print(t1-t0)
        cat("--> PE diff:", PE_tah, ", and PM diff:", PM_tah, ", and Sig diff:", Sig_tah, "\n\n")
        
        print(mu)
        cat("\n")
        print(table(uni_cell_types))
        cat("\n")
    }
    

    stopCluster(cl)
    
    PME_on_cell = t( apply(PM_on_cell*PE_on_cell,1,function(x)x/sum(x)) )
    PME_uni_cell = apply(PME_on_cell, 2, function(x) tapply(x, cell_id, max) )
    PME_uni_cell = PME_uni_cell[ match( cell_id, rownames(PME_uni_cell) ) , ]
    cell_types = colnames(PME_uni_cell)[apply(PME_uni_cell, 1, which.max)]
    names(cell_types) = rownames(PME_uni_cell)
    
    uni_cell_types = tapply( cell_types, names(cell_types), function(x) x[1] )
    count_cell_types = table(uni_cell_types)
    
    #cat( "   calculate adjacent similarity ... \n")
    #adjacent_similarity = get_adjacent_similarity(cells_on_spot, features, cell_types) 
    
    result <- list(lambda=lambda,
         mu=mu, 
         sigma=sigma, 
         PE_on_spot = PE_on_spot,
         PM_on_cell = PM_on_cell,
         PME_uni_cell = PME_uni_cell, 
         cell_types=cell_types,
         uni_cell_types=uni_cell_types,
         Signature=Signature,
         cells_on_spot=cells_on_spot)
    
    return(result)
}


