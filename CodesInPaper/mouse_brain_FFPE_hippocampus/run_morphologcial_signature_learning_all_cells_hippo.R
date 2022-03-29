deconvolution = TRUE
source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_hippocampus/parameters_hippo.R")

dwls = function(S, B)
{
    library(Matrix)
    source("/archive/SCCC/Hoshida_lab/s184554/Project/stRNAseq/DWLS/Deconvolution_functions.R")
    
    nb = ncol(B)
    ns = ncol(S)
    prop_mat <- as.data.frame(matrix(0, nrow = nb, ncol = ns) )
    rownames(prop_mat) <- colnames(B)
    colnames(prop_mat) <- colnames(S)
    for (ib in 1:nb) {
        print( sprintf("Estimating proportion for spot : %d / %d", ib, nb) )
        b <- B[,ib]
        tr <- trimData(S, b)
        tr$sig <- tr$sig[,colSums(tr$sig) > 0]
        is_pd <- eigen(t(tr$sig)%*%tr$sig)$values
        is_pd <- all(is_pd > 10e-6)
        if (!(is_pd)) { next}
        try(solDWLS <- solveDampenedWLS(tr$sig,tr$bulk),next)
        print("Proportions >> ")
        prop_mat[ib,names(solDWLS)] <- solDWLS
    }
    
    t(prop_mat)
    
}

############################################################
## deconvolution on multiorgan
############################################################

cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.8*spot_radius)

result <- STIE(ST_expr, Signature, cells_on_spot, features, 
               lambda=0, steps=30, 
               known_signature=TRUE, known_cell_types=FALSE)

################################################################################

cells_on_spot$cell_types = result$cell_types
spot_id = as.character(cells_on_spot$spot)
cell_id = as.character(cells_on_spot$cell_id)
Signature = Signature[rownames(Signature)%in%colnames(ST_expr),]
ST_expr2 = t( ST_expr[ rownames(ST_expr)%in%spot_id, match(rownames(Signature),colnames(ST_expr)) ] )

coefs = table( as.character(cells_on_spot$spot), as.character(cells_on_spot$cell_types) )
coefs = coefs[match( colnames(ST_expr2), rownames(coefs) ), ]
Signature2 = t( apply( ST_expr2, 1, function(x) solveOLS( coefs, as.matrix(x), scaled=F ) ) )

cc = data.frame( celltypes = c( "CA1", "CA2",  "CA3",    "DG",      "Glia"), 
                      cols = c( "red", "blue", "green", "black",  "yellow") )

Signature = Signature[, intersect(as.character(cc$celltypes),colnames(Signature)) ]
Signature2 = Signature2[, intersect(as.character(cc$celltypes),colnames(Signature2)) ]

prop2 = apply( Signature2, 2, function(x) solveNNLS(Signature,x) )  
corr2 = cor(log2(Signature2+1),log2(Signature+1))
prop_dwls2 = dwls(S=Signature, B=Signature2)

par(mfrow=c(2,2))
barplot(prop2, col=as.character(cc$cols), las=3, border='white')
barplot(prop_dwls2, col=as.character(cc$cols), las=3, border='white')
heatmap(corr2,Colv=NA,Rowv=NA)




