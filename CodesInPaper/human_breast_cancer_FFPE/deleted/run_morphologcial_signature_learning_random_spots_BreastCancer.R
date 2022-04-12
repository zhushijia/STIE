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

STIE_result = result

cells_on_spot = STIE_result[['2.8']]$cells_on_spot
cell_types = STIE_result[['2.8']]$cell_types

if(0)
{
    Signature = STIE_result[['2.8']]$Signature
    ST_expr = ST_expr[, match( rownames(Signature), colnames(ST_expr) )]
}

lambda=0
steps=20
known_signature=FALSE
known_cell_types=TRUE
min_cells=2
morphology_steps=10


############################################################
## clustering
############################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
load("HumanBreastCancer_clustering.RData")

STIE_result = result

cells_on_spot = STIE_result[[6]]$cells_on_spot
cell_types = STIE_result[[6]]$cell_types
lambda=0
steps=30
known_signature=FALSE
known_cell_types=TRUE
min_cells=2
morphology_steps=10

################################################################################
####### randomly assign cell types
################################################################################

num = 10

set.seed(1234)
selected_spots = sample( unique(as.character(cells_on_spot$spot)), num )
cells_on_spot$cell_types = NA
index = which( as.character(cells_on_spot$spot) %in% selected_spots )
cells_on_spot$cell_types[ index ] = cell_types[ index ]
table(cell_types[ index ])

coefs = table( as.character(cells_on_spot$spot)[ index ], as.character(cells_on_spot$cell_types)[ index ] )
ST_expr2 = t( ST_expr[ match( rownames(coefs), rownames(ST_expr)) ,] )
Signature = t( apply( ST_expr2, 1, function(x) solveOLS( coefs, as.matrix(x), scaled=F ) ) )

cor(STIE_result[[6]]$Signature, Signature)

Signature_ = Signature
ST_expr2 = t( ST_expr[ rownames(ST_expr)%in%spot_id, match(rownames(Signature),colnames(ST_expr)) ] )

mu = do.call(cbind, lapply(features, function(f) tapply( cells_on_spot[,f], 
                                                         as.character(cells_on_spot[,"cell_types"]), mean ) ))
sigma = do.call(cbind, lapply(features, function(f) tapply( cells_on_spot[,f], 
                                                            as.character(cells_on_spot[,"cell_types"]), sd ) ))
colnames(mu) = colnames(sigma) = features

PE_on_spot = t( apply( ST_expr2, 2, function(x) solveOLS( Signature, as.matrix(x) ) ) )
PM_on_cell = calculate_morphology_probability(cells_on_spot, features, mu, sigma )


################################################################################

result = STIE(ST_expr, Signature=Signature, cells_on_spot, features, lambda=lambda, steps=steps, 
              morphology_steps=morphology_steps,
              known_signature=FALSE, known_cell_types=TRUE)

table(STIE_result[['2.8']]$cell_types, result$cell_types)
table(STIE_result[[6]]$cell_types, result$cell_types)


