deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

calculate_mse = function( B, X=Signature, Y=t(ST_expr) ) {
    Y = Y[, colnames(Y) %in% rownames(B)]
    Y = Y[, match(rownames(B), colnames(Y))]
    Y = Y[match(rownames(X), rownames(Y)), ]
    X = X[, colnames(X) %in% colnames(B)]
    B = B[, match(colnames(X), colnames(B))]
    mse = sapply(1:ncol(Y), function(i) {
        y = as.numeric(Y[, i])
        b = as.numeric(B[i, ])
        bx = X %*% b
        t = sum(y)/sum(bx)
        mean((y - t * bx)^2)
    })
    mse
}

########################################################################
######## STIE
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
load("BreastCancer_spot_BIC_new_get_cells_on_spot.RData")
score = lapply(result, function(x) calculate_BIC(x,ST_expr)) 
stie_rmse = sapply(score,function(x) mean(sqrt(x$mse)))

########################################################################
######## DWLS
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/Visium_FFPE_Human_Breast_Cancer/DWLS_on_Wu_etal_2021_BRCA_Signature")
x = load("Visium_FFPE_Human_Breast_Cancer_DWLS_on_Wu_etal_2021_BRCA_Signature.RData" )
prop_mat = prop_mat[, !colnames(prop_mat)%in%c("Myeloid","PVL")]
dwls_mse = calculate_mse( B=prop_mat, X=Signature, Y=t(ST_expr) )
dwls_rmse = mean(sqrt(dwls_mse))
    
########################################################################
######## SPOTlight
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/SPOTlight")
spotlight_ls <- readRDS(file = "spotlight_ls.rds" )
nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]
rownames(decon_mtrx) = rownames(ST_expr)

spotlight_mse = calculate_mse( B=decon_mtrx[,-ncol(decon_mtrx)], X=Signature, Y=t(ST_expr) )
spotlight_rmse = mean(sqrt(spotlight_mse))



########################################################################
mse = lapply(score,function(x)x$mse)
mse[[length(mse)+1]] = dwls_mse
mse[[length(mse)+1]] = spotlight_mse
boxplot(mse,outline=F)

rmse = c(stie_rmse, dwls_rmse, spotlight_rmse)
plot(rmse)


