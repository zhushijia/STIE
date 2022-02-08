#' get_cellchat
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
get_cellchat <- function(STIE_result, 
                         database=c("human","mouse"), 
                         db_category=c("Secreted Signaling"), 
                         DEG_pthres=1e-5)
{
    # https://github.com/sqjin/CellChat
    
    library(quadprog)
    library(CellChat)
    library(data.table)
    library(patchwork)
    options(stringsAsFactors = FALSE)
    Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
    
    
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
    
    if(database=='human')  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
    if(database=='mouse')  CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
    #showDatabaseCategory(CellChatDB)
    
    if(is.null(db_category)) {
        CellChatDB.use <- CellChatDB
    } else {
        CellChatDB.use <- subsetDB(CellChatDB, search=db_category) # use Secreted Signaling
    }
    
    if(0)
    {
        gene.use_input <- extractGene(CellChatDB.use)
        complex_input <- CellChatDB.use$complex
        interaction_input <- CellChatDB.use$interaction
        pairLR <- select(interaction_input, ligand, receptor)
        LRgenes = setdiff( unique(c( unlist(complex_input), unlist(pairLR))), "" )
    }
    
    ST_expr2 = t( ST_expr[ rownames(ST_expr)%in%spot_id, match(rownames(Signature),colnames(ST_expr)) ] )
    ST_expr3 = t( ST_expr[ rownames(ST_expr)%in%spot_id, ] )
    
    coefs = do.call(rbind, lapply( 1:nrow(PE_on_spot), function(i) {
        
        cat(i,"\n")
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
    
    Signature3 = t( apply( ST_expr3, 1, function(x) solveOLS( coefs, as.matrix(x), scaled=F ) ) )
    
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
        y = log2(ST_expr3[i,]+1)
        z = summary(lm( y ~ coefs-1 ))$coef
        m0 = mean(y)
        sd0 = sd(y)
        n = colSums(coefs)
        
        info = do.call(rbind, lapply( 1:ncol(Signature), function(j) {
            with( ttest.2sample( m1=z[j,1], sd1=sqrt( z[j,2]^2*n[j] ), n1=n[j], 
                                m2=m0, sd2=sd0, n2=length(y) ), 
                 data.frame( clusters = colnames(Signature)[j],
                             features = rownames(ST_expr3)[i],
                             logFC = m,
                             pct.1 = m1,
                             pct.2 = m2,
                             pvalues = p))
            
        }))
        
        subset(info, pvalues<DEG_pthres)
        
    } ))
    
    features = unique(as.character(features.info$features))
    var.features = list(features.info=features.info,
                        features = features)
    
    data.input = Signature3[features,]
    meta = data.frame(cluster=colnames(Signature3),labels=colnames(Signature3))
    data.input = do.call(cbind, lapply(1:2, function(i) data.input ) )
    meta = do.call(rbind, lapply(1:2, function(i) meta ) )
    
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
    
    cellchat <- addMeta(cellchat, meta = meta)
    cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
    
    cellchat@DB <- CellChatDB.use
    
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    #future::plan("multiprocess", workers = 4) # do parallel
    
    cellchat@var.features = var.features
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    if(database=='human')  cellchat <- projectData(cellchat, PPI.human)
    if(database=='mouse')  cellchat <- projectData(cellchat, PPI.mouse)
    
    cellchat <- computeCommunProb(cellchat)
    #cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    groupSize <- as.numeric(table(cellchat@idents))
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    
    cellchat
}

