
############################################################
######## 
############################################################

layers_res_s2_0 = STIE(ST_expr, Signature, cell_coordinates=morphology_fts, features, spot_coordinates, spot_radius=2*spot_radius, lambda=0)
layers_res_s2_1e6 = STIE(ST_expr, Signature, cell_coordinates=morphology_fts, features, spot_coordinates, spot_radius=2*spot_radius, lambda=1e6)

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/parameter_testing")
save(layers_res_s2_0, file="layers_res_s2_0.RData")
parameter_testing(result_name="layers_res_s2_0", im, w=3000, h=3000, margin=100, spot_coordinates)


largeGroup_res_s3_0 = STIE(ST_expr, Signature, cell_coordinates=morphology_fts, features, spot_coordinates, spot_radius=3*spot_radius, lambda=0)
largeGroup_res_s2_1e6 = STIE(ST_expr, Signature, cell_coordinates=morphology_fts, features, spot_coordinates, spot_radius=2*spot_radius, lambda=1e6)
save(largeGroup_res_s3_0, file="largeGroup_res_s3_0.RData")
save(largeGroup_res_s2_1e6, file="largeGroup_res_s2_1e6.RData")
parameter_testing(result_name="largeGroup_res_s3_0", im, w=3000, h=3000, margin=100, spot_coordinates)
parameter_testing(result_name="largeGroup_res_s2_1e6", im, w=3000, h=3000, margin=100, spot_coordinates)





############################################################
# lambda
############################################################
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2*spot_radius)

res = list()
las = c(0,1,10,1e2,1e3,1e4,1e5,1e6)
for(i in 1:length(las))
{
    cat(las[i],'\n')
    res[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=las[i], steps=30, 
                    morphology_steps=ceiling(steps/3), 
                    known_signature=TRUE, known_cell_types=FALSE)
}

score1 = lapply( res, function(x) BIC(x) )
score2 = sapply(score1,mean)
plot(score2)

############################################################
# spot
############################################################
res2 = list()
rr = c(0.5,seq(1,3,0.2))
for(i in 1:length(rr))
{
    cat(rr[i],'\n')
    cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, rr[i]*spot_radius)
    res2[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=0, steps=30, 
                     morphology_steps=ceiling(steps/3), 
                     known_signature=TRUE, known_cell_types=FALSE)
}

score1 = lapply( res2, function(x) BIC(x) )
score2 = sapply(score1,mean)
names(score2) = names(score2) = rr
plot(score2)
rr[which.min(score2)]

############################################################
# cluster on all genes
############################################################
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2*spot_radius)

res3 = list()
for(i in 2:10)
{
    cat(i,'\n')
    
    cluster = read.csv(paste0("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/Visium_FFPE_Mouse_Brain/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
    all( as.character(cluster[,1]) %in% rownames(ST_expr) )
    ST_expr3 = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
    Signature = t(apply(ST_expr3, 2, function(x) tapply(x,cluster[,2],mean) ))
    
    res3[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=0, steps=30, 
                     morphology_steps=ceiling(steps/3), 
                     known_signature=FALSE, known_cell_types=FALSE)
}

score1 = lapply( res3[2:9], function(x) BIC(x) )
score2 = sapply(score1,mean)
names(score2) = names(score2) = rr
plot(score2)
rr[which.min(score2)]


############################################################
# cluster on known signature genes
############################################################
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2*spot_radius)

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/BroadInstitute_SingleCell")
Signature = read.delim("Major_cell_types_marker_genes.txt", header=T, row.names=1)
SS = Signature[,!colnames(Signature)%in%c("Ependymal.cells")]
SS = as.matrix(SS)[rownames(SS)%in%colnames(ST_expr),]

res4 = list()
for(i in 2:10)
{
    cat(i,'\n')
    
    
    cluster = read.csv(paste0("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/Visium_FFPE_Mouse_Brain/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
    all( as.character(cluster[,1]) %in% rownames(ST_expr) )
    ST_expr3 = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
    Signature = t(apply(ST_expr3, 2, function(x) tapply(x,cluster[,2],mean) ))
    Signature = Signature[rownames(Signature)%in%rownames(SS), ]
    
    res4[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=0, steps=30, 
                     morphology_steps=ceiling(steps/3), 
                     known_signature=FALSE, known_cell_types=FALSE)
}

score1 = lapply( res4[2:10], function(x) BIC(x) )
score2 = sapply(score1,mean)
names(score2) = names(score2) = rr
plot(score2)



setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/parameter_testing")
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/parameter_testing")
save(test, file="test.RData")
parameter_testing(result_name="test", im, w=3000, h=3000, margin=100, spot_coordinates)



