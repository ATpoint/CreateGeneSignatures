args <- commandArgs(trailingOnly=TRUE)
run_test <- if(is.na(args[1])) FALSE else as.logical(args[1])

if(run_test){
  
  if(!requireNamespace(c("edgeR", "pheatmap")))
    BiocManager::install("edgeR", "pheatmap")
  
  library(edgeR)
  library(pheatmap)

  source("R/RankDEGs.R")
  source("R/CreateGeneSignatures.R")
  
  res <- readRDS("inst/extdata/res.rds")
  
  #/ rank DEGs
  ranked <- RankDEGs(res)
  
  #/ test1: standard use
  signatures <- CreateGeneSignatures(ranked=ranked, keep.n=50, min.prop=1, extended=FALSE)
  
  #/ test2: pool groups
  signatures2 <- CreateGeneSignatures(ranked=ranked, use_groups=c("CD4T", "CD8T"), extended=FALSE)
  
  #/ test3: exclude groups
  signatures3 <- CreateGeneSignatures(ranked=ranked, exclude_groups=c("CD4T"), extended=FALSE)
  
  #/ test 4: with extended
  signatures4 <- CreateGeneSignatures(ranked=ranked, min.prop=2/3, extended=TRUE)

  
  
  counts <- readRDS(paste0(
    system.file("extdata",package="CreateGeneSignatures"),
    "/haemopedia_subset.rds"))
  
  y <- DGEList(counts=counts,group=gsub("\\..", "", colnames(counts)))
  
  design <- model.matrix(~0+group,y$samples)
  
  colnames(design) <- gsub("group", "", colnames(design))
  
  y <- y[filterByExpr(y),,keep.lib.size=FALSE]
  
  y <- calcNormFactors(y)
  
  y <- estimateDisp(y,design)
  
  fit <- glmQLFit(y,design)
  
  contrasts <- makeContrasts(CD4T_vs_CD8T  = CD4T-CD8T,
                             CD4T_vs_NK    = CD4T-NK,
                             CD4T_vs_NveB  = CD4T-NveB,
                             CD8T_vs_NK    = CD8T-NK,
                             CD8T_vs_NveB  = CD8T-NveB,
                             NK_vs_NveB    = NK-NveB,
                             levels = design)
  
  res <- sapply(colnames(contrasts), function(con){
    tt<-topTags(glmTreat(fit,contrast=contrasts[,con],lfc=log2(1.5)),n=Inf)$table
    return(data.frame(Gene=rownames(tt), tt)[,c("Gene", "logFC", "PValue", "FDR")])
  }, simplify = FALSE)
  
  logcpm <- log2(edgeR::cpm(y,log=FALSE)+1)
  
  col_order <- unlist(lapply(names(ranked), 
                             function(x) grep(paste0("^", x), colnames(logcpm))))
  
  logcpmZ <- t(scale(t(logcpm[unique(unlist(signatures)),])))
  
  p <- pheatmap(mat=logcpmZ[,col_order],
                show_rownames=FALSE, cluster_rows=FALSE, cluster_cols=FALSE)
  
  png("heatmap.png")
  print(p); dev.off()

}
