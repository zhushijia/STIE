simulationData = function(sd_mr, sd_er)
{    
    library(STIE)
    library(quadprog)
    
    set.seed(1234)
    
    spot_num = 1000
    gene_num = 500
    celltype_num = 10
    cell_num = 10000
    sc_num = 100
    feature_num = 3
    mu = matrix( runif( celltype_num*feature_num, 0, 3 ), nrow=celltype_num, ncol=feature_num )
    sigma = matrix( runif( celltype_num*feature_num, 0, sqrt(sd_mr) ), nrow=celltype_num, ncol=feature_num )
    rownames(mu) = rownames(sigma) = paste0("type",1:celltype_num)
    colnames(mu) = colnames(sigma) = paste0("feature",1:feature_num)
    
    features = paste0("feature",1:feature_num)
    cell_id = paste0( 'cell', 1:cell_num )
    spot_id = sample( paste0("spot",1:spot_num), cell_num, replace = TRUE )
    cell_types = sample( paste0("type",1:celltype_num), cell_num, replace = TRUE )
    
    Signature = matrix(data=0, nrow=gene_num, ncol=celltype_num)
    colnames(Signature) = paste0("type",1:celltype_num)
    rownames(Signature) = paste0("gene",1:gene_num)
    
    morphology_fts = matrix(data=0, nrow=cell_num, ncol=feature_num)
    colnames(morphology_fts) = features
    cells_on_spot <- data.frame(cell_id=cell_id, spot=spot_id, cell_types=cell_types )
    
    ST_expr =   matrix(data=0, nrow=spot_num, ncol=gene_num)
    rownames(ST_expr) = paste0("spot",1:spot_num)
    colnames(ST_expr) = paste0("gene",1:gene_num)
    
    cell_prop = table(spot_id, cell_types)
    cell_prop = cell_prop[ match( rownames(ST_expr), rownames(cell_prop) ), 
                           match( colnames(Signature), colnames(cell_prop) )]
    
    for(i in 1:ncol(Signature)) Signature[,i] = runif(gene_num,0,celltype_num)  #rnbinom(gene_num, mu=i, size=1)
    ST_expr = cell_prop%*%t(Signature) 
    for(i in 1:nrow(ST_expr)) 
        for(j in 1:ncol(ST_expr))
            ST_expr[i,j] = rnbinom(1, mu=ST_expr[i,j], size=1/(2*sd_er) )
    
    
    for(i in 1:nrow(mu)) {
        indi = which( cell_types==rownames(mu)[i] )
        morphology_fts[indi, ] = do.call(cbind, lapply(1:ncol(mu), function(j) {
            rnorm(length(indi),mu[i,j],sigma[i,j])
        }))
    }
    
    cells_on_spot = data.frame(cells_on_spot, morphology_fts)
    

    SC_expr = do.call(cbind, lapply( 1:celltype_num, function(i) {
      x = Signature[,i]
      do.call(cbind,lapply(1:sc_num, function(j) x ))
    }))
    colnames(SC_expr) = paste0("cell",1:(celltype_num*sc_num) )
    rownames(SC_expr) = paste0("gene",1:gene_num)
    for(i in 1:nrow(SC_expr)) 
      for(j in 1:ncol(SC_expr))
        SC_expr[i,j] = rnbinom(1, mu=SC_expr[i,j], size=1/1)
    
    annot = data.frame( cell_id=colnames(SC_expr), cell_type=rep(colnames(Signature),each=sc_num) )
    rownames(annot) = colnames(SC_expr)
    
    ##############################################################################
    ###### only use gene expression and morphological features to estimate
    ##############################################################################
    pred_prop_expr = t(apply(ST_expr,1,function(x)solveNNLS(S=Signature, B=x, scaled=T)))
    
    km = kmeans(morphology_fts,celltype_num)$cluster
    km = kmeans(morphology_fts,celltype_num)$cluster
    pred_prop_morph = table(cells_on_spot$spot,km)
    pred_prop_morph = pred_prop_morph[rownames(cell_prop),]
    pred_prop_morph = t(apply(pred_prop_morph,1,function(x)x/sum(x)))
    
    ind = c()
    for(i in 1:celltype_num) {
      cc = cor(pred_prop_morph, cell_prop[,i])
      ind[i] = setdiff(order(cc,decreasing=T),ind)[1]
    }
    
    pred_prop_morph = pred_prop_morph[,ind]
    
    #total_cc = sapply(1:celltype_num,function(i) cor(pred_prop_expr[,i],pred_prop_morph[,i]) )
    corr_expr = sapply(1:nrow(pred_prop_expr),function(i)cor(cell_prop[i,],pred_prop_expr[i,]))
    corr_morph = sapply(1:nrow(pred_prop_morph),function(i)cor(cell_prop[i,],pred_prop_morph[i,]))
    
    corr_inter = sapply(1:nrow(pred_prop_expr),function(i) cor(as.numeric(pred_prop_expr[i,]),as.numeric(pred_prop_morph[i,])) )
    
    corr = list(corr_expr=corr_expr,corr_morph=corr_morph,corr_inter=corr_inter)
    boxplot(corr,las=3,ylim=c(0,1))
    
    cat("sd_mr=", sd_mr, "sd_er=", sd_er, "\n")
    
    ##############################################################################
    mse_expr = sapply(1:nrow(ST_expr), function(i) {
        y = as.matrix( ST_expr[i,] )
        b = cell_prop[i,]
        bx = Signature%*%b
        t = sum(y)/sum(bx)
        r = (y-t*bx)
        mean(r^2)
    } )
    
    ##############################################################################
    mu_ = apply(morphology_fts,2,function(x)tapply(x,cell_types,mean))
    mu_ = mu_[rownames(mu), colnames(mu)]
    sigma_ = apply(morphology_fts,2,function(x)tapply(x,cell_types,sd))
    sigma_ = sigma_[rownames(mu), colnames(mu)]
    
    PM = matrix(nrow = nrow(morphology_fts), ncol = nrow(mu_), data = 0)
    rownames(PM) = cell_id
    colnames(PM) = rownames(mu_)
    for (t in 1:nrow(mu_)) {
        PM[, t] <- apply(as.matrix(morphology_fts), 1, function(x) 
            prod(dnorm(as.numeric(x), mean = as.numeric(mu_[t,]), sd = as.numeric(sigma_[t, ]))))
    }
    PM[PM < 0] = 1e-100
    PME_on_cell_unscaled = PM*(1/ncol(PM))
    ll_morph = log(rowSums(PME_on_cell_unscaled))/ncol(mu)
    
    ##############################################################################
    fileName = paste0("SimulationData_sd_mr_",sd_mr, "_sd_er_", sd_er, ".RData")
    save( ST_expr, cells_on_spot, morphology_fts, features, 
          mu, sigma, 
          Signature, SC_expr, annot,
          cell_types, cell_prop, 
          pred_prop_expr, pred_prop_morph, 
          corr_expr, corr_morph, corr_inter,
          mse_expr=mse_expr, ll_morph=ll_morph,
          file=fileName )
    ##############################################################################
    
}

simulationSummary = function(result)
{
    x= table(cells_on_spot$cell_types, result$cell_types)
    
    print( sum(diag(x)/sum(x)) )
    print( cor( mu[rownames(result$mu),], result$mu ) )
    print( cor( sigma[rownames(result$mu),], result$sigma ) )
    
    pred_prop = table(cells_on_spot$spot, result$cell_types)
    pred_prop = pred_prop[ match( rownames(ST_expr), rownames(pred_prop) ), 
                           match( colnames(Signature), colnames(pred_prop) )]
    pred_prop = t(apply(pred_prop,1,function(x) x/sum(x,na.rm=T)))
    
    corr = sapply(1:nrow(pred_prop),function(i)
        cor.test(cell_prop[i,],pred_prop[i,],na.rm=T)$estimate)
    
    mse = sapply(1:nrow(ST_expr), function(i) {
        y = as.matrix( ST_expr[i,] )
        b = pred_prop[i,]
        bx = Signature%*%b
        t = sum(y)/sum(bx)
        r = (y-t*bx)
        mean(r^2)
    } )
    
    list(consistency=x, pred_prop=pred_prop, prop_corr=corr, mse=mse)
}

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/Simulation/data")
simulationData( sd_mr = 3, sd_er = 1 )
simulationData( sd_mr = 3, sd_er = 3 )
simulationData( sd_mr = 1, sd_er = 3 )
simulationData( sd_mr = 1, sd_er = 1 )

##############################################################################
##############################################################################
##############################################################################
# sd_mr = 3; sd_er = 1
# sd_mr = 3; sd_er = 3
# sd_mr = 1; sd_er = 3
# sd_mr = 1; sd_er = 1

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/Simulation/data")
fileName = paste0("SimulationData_sd_mr_",sd_mr, "_sd_er_", sd_er, ".RData")
x = load(fileName)

library(STIE)
library(quadprog)

las = c(0,10^c(1:6))
STIE_results <- list()
for( i in 1:length(las) )
{
    la = las[i]
    cat('lambda=', la, "\n")
    STIE_results[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, 
                        lambda=la, steps=30, 
                        known_signature=TRUE, known_cell_types=FALSE, 
                        min_cells=-1, equal_prior=F) 
}
names(STIE_results) = las
STIE_sums = lapply( STIE_results, simulationSummary )
STIE_corrs = lapply(STIE_sums,function(x) x$prop_corr )
boxplot( STIE_corrs, ylim=c(0-1,1), outline=F )
##############################################################################
parentDir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/Simulation/results"
workdir=paste0(parentDir,"/SimulationData_sd_mr_",sd_mr, "_sd_er_", sd_er)
dir.create(workdir)
setwd(workdir)
save(STIE_results, STIE_corrs, STIE_sums, file="STIE_result2.RData")
##############################################################################
sd_mr
sd_er


