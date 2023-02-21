simulationData = function(sd_mr, sd_er, gene_num, celltype_num, spot_num, cell_num)
{    
    library(STIE)
    library(quadprog)
    
    set.seed(1234567)
    
    #spot_num = 1000
    #gene_num = 500
    #celltype_num = 10
    #cell_num = 10000
    
    feature_num = 3
    sc_num = 100
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
    
    
    list(ST_expr=ST_expr, Signature=Signature, cells_on_spot=cells_on_spot, features=features)
}

spot_nums = c(100,200,500,1000,2000,5000,10000)
time_spot = list()
for(i in 1:length(spot_nums))
{   
    cat(spot_nums[i],'\n')
    simD = simulationData( sd_mr = 1, sd_er = 1, gene_num=1000, celltype_num=10, spot_num=spot_nums[i], cell_num=10*spot_nums[i])
    attach(simD)
    time_spot[[i]] = system.time( res <- STIE(ST_expr, Signature, cells_on_spot, features, 
                                       lambda=1e4, steps=30, 
                                       known_signature=TRUE, known_cell_types=FALSE, 
                                       min_cells=-1, equal_prior=F) )
    print(time_spot[[i]])
}



celltype_nums = c(3,5,10,15,20,25,30)
time_celltype = list()
for(i in 1:length(celltype_nums))
{   
    cat(celltype_nums[i],'\n')
    simD = simulationData( sd_mr = 1, sd_er = 1, gene_num=1000, celltype_num=celltype_nums[i], spot_num=1000, cell_num=10000)
    attach(simD)
    time_celltype[[i]] = system.time( res <- STIE(ST_expr, Signature, cells_on_spot, features, 
                                              lambda=1e4, steps=30, 
                                              known_signature=TRUE, known_cell_types=FALSE, 
                                              min_cells=-1, equal_prior=F) )
    print(time_celltype[[i]])
}


gene_nums = c(100,200,500,1000,2000,5000,10000)
time_gene = list()
for(i in 1:length(gene_nums))
{   
    cat(gene_nums[i],'\n')
    simD = simulationData( sd_mr = 1, sd_er = 1, gene_num=gene_nums[i], celltype_num=10, spot_num=1000, cell_num=10000)
    attach(simD)
    time_gene[[i]] = system.time( res <- STIE(ST_expr, Signature, cells_on_spot, features, 
                                                  lambda=1e4, steps=30, 
                                                  known_signature=TRUE, known_cell_types=FALSE, 
                                                  min_cells=-1, equal_prior=F) )
    print(time_gene[[i]])
}


parentDir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/Simulation/results"
setwd(parentDir)
save( spot_nums, time_spot, celltype_nums, time_celltype, gene_nums, time_gene, file='runningTime_onSimulation.RData' )

pdf("runningTime_onSimulation.pdf")
par(mfrow=c(2,2))
barplot(sapply(time_spot,function(x)x[3]), names=spot_nums, las=3, ylab='second',xlab='spot #',border=NA)
barplot(sapply(time_gene,function(x)x[3]), names=gene_nums, las=3, ylab='second',xlab='gene #',border=NA)
barplot(sapply(time_celltype,function(x)x[3]), names=celltype_nums, las=3, ylab='second',xlab='celltype #',border=NA)
dev.off()
