#' Create gene signatures in a rank-based style
#' 
#' Convert lists of ranked genes into celltype-specific signatures
#' 
#' @param ranked nested lists of ranked genes for each celltype, see details and the RankDEGs function
#' @param delim a string indicating the delimiter of names(ranked), e.g. celltype1_vs_celltype2 would be "_vs_"
#' @param keep.n number of genes to keep per signature, by default all candidates are returned for manual posthoc filtering
#' @param min.prop minimum proportion of comparisons per celltype that a gene must be included in. See details.
#' 
#' @author Alexander Toenges
#' 
#' @examples 
#' # load RNA-seq data for CD4T-, CD8T and naive B cells from Haemopedia:
#' counts <- readRDS(paste0(
#'             system.file("extdata",package="CreateGeneSignatures"),
#'             "/haemopedia_subset.rds"))
#'
#' # Use edgeR to perform all pairwise comparisons
#' 
#' library(edgeR)
#' y <- DGEList(counts=counts,group=gsub("\\..", "", colnames(counts)))
#' design <- model.matrix(~0+group,y$samples)
#' colnames(design) <- gsub("group", "", colnames(design))
#' y <- y[filterByExpr(y),,keep.lib.size=FALSE]
#' y <- calcNormFactors(y)
#' y <- estimateDisp(y,design)
#' fit <- glmQLFit(y,design)
#' 
#' # all unique pairwise contrasts:
#' contrasts <- makeContrasts(CD4T_vs_CD8T  = CD4T-CD8T,
#'                            CD4T_vs_NK    = CD4T-NK,
#'                            CD4T_vs_NveB  = CD4T-NveB,
#'                            CD8T_vs_NK    = CD8T-NK,
#'                            CD8T_vs_NveB  = CD8T-NveB,
#'                            NK_vs_NveB    = NK-NveB,
#'                            levels = design)
#'                            
#' # test using glmTreat against a minumum fold change of ~1.5                          
#' res <- sapply(colnames(contrasts), function(con){
#'   tt<-topTags(glmTreat(fit,contrast=contrasts[,con],lfc=log2(1.5)),n=Inf)$table
#' return(data.frame(Gene=rownames(tt), tt))
#' }, simplify = FALSE)
#' 
#' # Rank the DEGs:
#' ranked <- RankDEGs(res)
#' 
#' # Create signatures, keeping top 50 signature genes that separate the respective celltype
#' # from all other celltypes:
#' signatures <- CreateGeneSignatures(ranked=ranked, keep.n=50, min.prop=1)
#' # check number of genes. for CD8T cells we found < 50 genes:
#' lapply(signatures,length) 
#' 
#' # Inspect signatures using heatmaps plotting the scaled logcpms of the signature genes
#' library(pheatmap)
#' logcpm <- log2(edgeR::cpm(y,log=FALSE)+1)
#' 
#' # either as a whole, deactivating row clustering to preserve order from signatures
#' # but enable columns clustering:
#' pheatmap(mat=t(scale(t(logcpm[unique(unlist(signatures)),]))), 
#' show_rownames=FALSE, cluster_rows = FALSE)
#' 
#' # or each signature individually:
#' pheatmap(mat=t(scale(t(logcpm[signatures$CD4T,]))), show_rownames=FALSE)
#' pheatmap(mat=t(scale(t(logcpm[signatures$CD8T,]))), show_rownames=FALSE)
#' pheatmap(mat=t(scale(t(logcpm[signatures$NveB,]))), show_rownames=FALSE)
#' pheatmap(mat=t(scale(t(logcpm[signatures$NK,]))), show_rownames=FALSE)
#' 
#' The signatures represent those combination of genes separate each of the
#' individual celltypes from all other celltypes. Here the parameters were very strict,
#' e.g. min.prop=1. This might make sense if celltypes are very different from each other,
#' so there are in case sufficient numbers of genes per celltype being overexpressed compared
#' to all other celltypes. If cells are expected to be less distinct, e.g. hematopoietic
#' progenitors that are developmentally "close" e.g. in a developmental continuum/trajectory
#' one might need to lower min.prop to get a sufficient number of genes per signature.
#' The user is encouraged to try different settings, checking how many genes per celltype
#' can be obtained with the given setting followed by inspection of the separation using heatmaps
#' as described in the examples. In might be the case though that cells cannot be reliably distinguished
#' using transcriptomics alone, or that only very few or no signature genes can be obtained.
#' The latter could be the case if focusing on a celltype that is developmentally "between"
#' two cell types, say a celltype that comes developmentally after a stem cell and upstream of
#' a more differentiated progenitor. This celltype might express pretty much all genes of both aforementioned
#' celltypes but at slightly different levels, not enough though to quality as markers.
#' Plotting heatmaps using the pooled marker genes of all celltypes should then show that these cells
#' (at modest levels) express genes that are markers in the other celltypes.
#' 
#' One should remember that a signature as a whole should separate a celltype from all
#' other celltypes. Therefore not every gene by itself must be upregulated int he given 
#' celltype against all other celltypes but the combination of all signature genes is the key
#' for a proper celltype separation.
#' 
#' @details 
#' The function takes as input a nested list that, for every celltype, contains the ranked genes from pairwise comparisons against every other
#' celltype in this analysis, see also the \code{RankDEGs} function of this package. For every celltype first a rank matrix is created that
#' stores the ranked genes from every pairwise comparison, assigning low scores to genes that ranked high in the input ranking and vice versa.
#' Genes to form the signatures are then selected by choosing those genes with minimal ranks across all comparisons.
#' In an iterative fashion genes that ranked high all comparisons are chosen first (and then ordered by median score), then genes ranked high in all but one, 
#' in all but two, all but three etc comparisons. The user can choose via \code{min.prop} the fraction of comparisons that a gene be highly ranked in.
#' For example a setup with four celltypes and \code{min.prop=1} would mean that only genes are included that ranked high in all comparisons.
#' If \code{min.prop=2/3} is chosen then signature genes must rank high in at least two out of three comparisons. The higher this value between (between zero and one)
#' the more stringent the filtering. If celltypes are decently separated it might be trivial to find large numbers of signature genes even with min.prop of 1.
#' In cases where cells are very similar to each other, e.g. cells being close in a developmental continuum, one will need to lower this threshold in order to get some
#' (or any) signature genes. The user is encouraged to try different values and check the set of marker genes using heatmaps to inspect the spearation between celltypes
#' as suggested in the examples section. 
#' 
#' @export
CreateGeneSignatures <- function(ranked, delim="_vs_", keep.n=Inf, 
                                 min.prop=0.75){
  
  if(is.null(names(ranked))) stop("ranked has no names")
  if(!min.prop <= 1 & min.prop > 0) stop("min.prop must be between > 0 and <= 1")
  
  nms <- names(ranked)
  
  s <- sapply(nms, function(x){
    
    genes.ranged <- ranked[[x]]
    genes.ranged[[x]] <- NULL # if the list contains the reference cell type itself
    
    #/ rank matrix, low values means high ranks and viceversa, 1 means not present in that pairwise list: 
    rmat <- rankmatrix(genes.ranged = genes.ranged)
    
    #/ Now minimize the ranks in an iterative fashion.
    #/ Start with genes present in all comparison, n-1, n-2 etc.
    #/ Function is simply the mean
    isna <- !is.na(rmat)
    from=ncol(rmat)
    to=round(min.prop*from)
    med <- function(y) median(y,na.rm=TRUE)
    
    final <- unlist(lapply(from:to, function(n){
      
      idx <- base::rowSums(isna)==n
      if(sum(idx)==0) return(NULL)
      mat <- apply(rmat[idx,,drop=FALSE], 1, med)
      mat[order(mat)]
      
    }))
    
    return(head(names(final), keep.n))
    
  },simplify = FALSE)
  
  return(s)
  
}

#' For the genes being present in genes.ranged per celltype we rank from 1 (best rank) to length 
#' of the current genes.ranged entry, and NA for genes not in that respective entry but in u:
rankmatrix <- function (genes.ranged) 
{
  u = unique(c(genes.ranged, recursive = TRUE))
  N = length(u)
  
  #/ matrix with only NAs, to be filled with the ranks:
  rmat = matrix(NA, nrow=length(u), ncol=length(genes.ranged), 
                dimnames=list(u, names(genes.ranged)))
  
  #/ fill matrix with ranks per gene and celltype. I
  #/ If gene is missing in a celltype (=was not signif vs other celltypes),
  #/ then leave at NA:
  for (i in names(genes.ranged)) {
    
    rmat[genes.ranged[[i]], i] = seq(1, length(genes.ranged[[i]]))
    
  }
  return(rmat)
}

