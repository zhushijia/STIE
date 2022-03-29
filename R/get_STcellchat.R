#' get_STcellchat
#'
#' @param STIE_result 
#' @param ST_expr 
#' @param database 
#' @param db_category 
#' @param max_reps 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
#' 
get_STcellchat <- function(STIE_result, ST_expr, 
                         database=c("human","mouse"), 
                         db_category=c("Secreted Signaling"), 
                         max_reps=NULL )
{
    # https://github.com/sqjin/CellChat
    # http://www.cellchat.org/
    
    if(0)
    {
        #STIE_result = result
        database="human"
        db_category=NULL
        max_reps=1000
    }
    
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
    
    if(database=='human')  CellChatDB <- CellChatDB.human
    if(database=='mouse')  CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
    #showDatabaseCategory(CellChatDB)
    
    if(is.null(db_category)) {
        CellChatDB.use <- CellChatDB
    } else {
        CellChatDB.use <- subsetDB(CellChatDB, search=db_category) # use Secreted Signaling
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
    
    Signature3 = t( apply( ST_expr3, 1, function(x) solveNNLS( coefs, as.matrix(x), scaled=F ) ) )
    
    ######################################################################################################
    ############## pseudo-data
    ######################################################################################################
    
    uni_cellid = unique(cell_id)
    uni_cellid_celltypes = cell_types[ match( uni_cellid, names(cell_types) ) ]
    num_celltypes = table(uni_cellid_celltypes)
    
    reps = num_celltypes[ match( colnames(Signature3), names(num_celltypes) ) ]
    
    data.input = do.call(cbind, lapply( 1:ncol(Signature3), function(i) {
            ni = ifelse( is.null(max_reps), reps[i], min(reps[i],max_reps) )
            xi = Signature3[,i]
            do.call(cbind, lapply(1:ni, function(j) xi ) )
        } ) )
    colnames(data.input) = NULL
    
    meta = data.frame(cluster=rep(colnames(Signature3),reps),
                      labels=rep(colnames(Signature3),reps) )
    
    ######################################################################################################
    
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
    cellchat <- addMeta(cellchat, meta = meta)
    cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
    
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    cellchat <- identifyOverExpressedInteractions(cellchat, features=rownames(Signature3))
    
    if(database=='human')  cellchat <- projectData(cellchat, PPI.human)
    if(database=='mouse')  cellchat <- projectData(cellchat, PPI.mouse)
    
    cellchat <- computeCommunProb(cellchat)
    #cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    
    par(mfrow = c(1,2), xpd=TRUE)
    groupSize <- as.numeric(table(cellchat@idents))
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    
    if(0)
    {
        ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
        ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
        ht1[1:20,] + ht2[1:20,]
    }
    
    cellchat
    
}

