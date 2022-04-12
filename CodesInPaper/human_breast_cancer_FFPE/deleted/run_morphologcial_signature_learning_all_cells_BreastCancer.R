deconvolution = FALSE
clustering = FALSE
signature_learning = TRUE

known_signature = FALSE
known_cell_types = TRUE

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

############################################################
## deconvolution
############################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
load("BreastCancer_spot_BIC_new_get_cells_on_spot.RData")

STIE_result = result[['2.8']]

cells_on_spot = STIE_result$cells_on_spot
cell_types = STIE_result$cell_types
Signature = STIE_result$Signature
ST_expr = ST_expr[, match( rownames(Signature), colnames(ST_expr) )]
lambda=0
known_signature=FALSE
known_cell_types=TRUE
min_cells=2
steps=20
morphology_steps=10



################################################################################
####### randomly assign cell types
################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq")
y = load("Wu_etal_2021_BRCA_scRNASeq_DWLS_Signature.RData")
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/Visium_FFPE_Human_Breast_Cancer/DWLS_on_Wu_etal_2021_BRCA_Signature")
x = load("Visium_FFPE_Human_Breast_Cancer_DWLS_on_Wu_etal_2021_BRCA_Signature.RData" )

################################################################################
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2*spot_radius)
spot_id = as.character(cells_on_spot$spot)
cell_id = as.character(cells_on_spot$cell_id)
Signature = Signature[rownames(Signature)%in%colnames(ST_expr),]
ST_expr2 = t( ST_expr[ rownames(ST_expr)%in%spot_id, match(rownames(Signature),colnames(ST_expr)) ] )

coefs = table( as.character(cells_on_spot$spot), as.character(cells_on_spot$cell_types) )
coefs = coefs[match( colnames(ST_expr2), rownames(coefs) ), ]
Signature2 = t( apply( ST_expr2, 1, function(x) solveOLS( coefs, as.matrix(x), scaled=F ) ) )
apply( Signature2, 2, function(x) solveOLS(Signature,x) )  



result = STIE(ST_expr, Signature=Signature, cells_on_spot, features, lambda=lambda, steps=steps, 
              morphology_steps=morphology_steps,
              known_signature=FALSE, known_cell_types=TRUE)

table(STIE_result$cell_types, result$cell_types)
cor(STIE_result$Signature, result$Signature, method='spearman')

apply( result$Signature, 2, function(x) solveOLS(STIE_result$Signature,x) )  
apply( Signature, 2, function(x) solveOLS(STIE_result$Signature,x) )  

################################################################################

if(0)
{
    im_scaled <- read_image(image=args$image, x_scale=args$x_scale)
    
    cells_on_spot$pixel_x = cells_on_spot$X * args$x_scale
    cells_on_spot$pixel_y = cells_on_spot$Y * args$x_scale
    uni_celltypes = sort(unique(as.character(cells_on_spot$cell_types)))
    cols = c("red", "purple", "green","cyan", "yellow", "black", "blue")
    
    #setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
    #pdf("Breast_Cancer_FFPE_DL_celltypes_information_model40_image20.pdf")
    display(im_scaled, method='raster')
    for(i in 1:length(uni_celltypes))
    {
        X = subset( cells_on_spot, cell_types==uni_celltypes[i] )[,c('pixel_x','pixel_y')]
        #X = subset( cells_on_spot, result$cell_types==uni_celltypes[i] )[,c('pixel_x','pixel_y')]
        points(X,cex=0.1,col=cols[i])
    }
    plot(NA,xlim=c(0,10),ylim=c(0,10))
    legend('topright',legend=uni_celltypes,col=cols,pch=1)
    
    #dev.off()
    
}


apply( result$Signature, 2, function(x) solveOLS(Signature,x) )  



