
source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

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

################################################################################
####### randomly assign cell types
################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq")
y = load("Wu_etal_2021_BRCA_scRNASeq_DWLS_Signature.RData")
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/Visium_FFPE_Human_Breast_Cancer/DWLS_on_Wu_etal_2021_BRCA_Signature")
x = load("Visium_FFPE_Human_Breast_Cancer_DWLS_on_Wu_etal_2021_BRCA_Signature.RData" )

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

cc = data.frame( celltypes = c( "Plasmablasts", "Bcells",  "Tcells",    "Myeloid",      "CAFs",  "Endothelial",  "PVL", "CancerEpithelial", "NormalEpithelial"), 
                      cols = c( "cyan",    "steelblue", "darkgreen", "darkorange",  "purple",  "darkred",     "yellow",   "black", "grey") )

Signature = Signature[, intersect(as.character(cc$celltypes),colnames(Signature)) ]
Signature2 = Signature2[, intersect(as.character(cc$celltypes),colnames(Signature2)) ]

prop2 = apply( Signature2, 2, function(x) solveOLS(Signature,x) )  
corr2 = cor(log2(Signature2+1),log2(Signature+1))
prop_dwls2 = dwls(S=Signature, B=Signature2)

barplot(prop2, col=as.character(cc$cols), las=3, border='white')
barplot(prop_dwls2, col=as.character(cc$cols), las=3, border='white')
heatmap(corr2,Colv=NA,Rowv=NA)

############################################################
## deconvolution on HD staining 
############################################################

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/")
cells = data.frame(fread("Breast_Cancer_FFPE_information_model40_image20.csv"))
morphology_fts = with(cells, data.frame(
    cell_id = paste0("cell_",rownames(cells)),
    cell_types=nucleus_class,
    probability=probability,
    X=centroid_x0,
    Y=centroid_y0,
    Area=area,
    Major=major_axis_length,
    Minor=minor_axis_length,
    Eccentricity=eccentricity,
    Orientation=orientation,
    Solidity=solidity) )
#morphology_fts = subset(morphology_fts, ! cell_types %in% c("blood") )
morphology_fts = morphology_fts[order(morphology_fts$probability,decreasing=T),]
morphology_fts = subset(morphology_fts,probability>0.5)
morphology_fts$pixel_x = morphology_fts$X * args$x_scale
morphology_fts$pixel_y = morphology_fts$Y * args$x_scale

############################################################
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.8*spot_radius)
spot_id = as.character(cells_on_spot$spot)
cell_id = as.character(cells_on_spot$cell_id)
Signature = Signature[rownames(Signature)%in%colnames(ST_expr),]
ST_expr2 = t( ST_expr[ rownames(ST_expr)%in%spot_id, match(rownames(Signature),colnames(ST_expr)) ] )

coefs = table( as.character(cells_on_spot$spot), as.character(cells_on_spot$cell_types) )
coefs = coefs[match( colnames(ST_expr2), rownames(coefs) ), ]
Signature3 = t( apply( ST_expr2, 1, function(x) solveNNLS( coefs, as.matrix(x), scaled=F ) ) )

Signature3 = Signature3[, c( "blood", "lymphocyte", "macrophage", "stroma", "tumor", "ductal epithelium", "necrosis" ) ]

prop3 = apply( Signature3, 2, function(x) solveOLS(Signature,x) )  
corr3 = cor(log2(Signature3+1),log2(Signature+1))
prop_dwls3 = dwls(S=Signature, B=Signature3)

par(mfrow=c(2,2))
barplot(prop3, col=as.character(cc$cols), las=3, border='white')
barplot(prop_dwls3, col=as.character(cc$cols), las=3, border='white')
heatmap(corr3,Colv=NA,Rowv=NA)





