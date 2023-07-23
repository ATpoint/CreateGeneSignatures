#' Create gene signatures in a rank-based style
#' 
#' Convert lists of ranked genes into celltype-specific signatures
#' 
#' @param ranked nested lists of ranked genes for each celltype, see details and the RankDEGs function
#' @param use_groups trigger group-aware mode. A vector with celltype names to be used to define markers for as a group,
#' so this group versus all others. See details.
#' @param exclude_groups vector, remove these groups from the signature creation
#' @param delim a string indicating the delimiter of names(ranked), e.g. celltype1_vs_celltype2 would be "_vs_"
#' @param keep.n number of genes to keep per signature, by default all candidates are returned for manual posthoc filtering
#' @param min.prop minimum proportion of comparisons per celltype that a gene must be included in. See details.
#' @param extended logical, whether to output the signature genes as a data.frame that indicates whether the gene
#' qualified as a marker against the other groups (1) or not (0)
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
#' signatures <- CreateGeneSignatures(ranked=ranked, keep.n=50)
#' # check number of genes. for CD8T cells we found < 50 genes:
#' lengths(signatures)
#' 
#' 
#' # Inspect signatures using heatmaps plotting the scaled logcpms of the signature genes
#' library(pheatmap)
#' logcpm <- log2(edgeR::cpm(y,log=FALSE)+1)
#' 
#' # plot a heatmap in the order of names(signatures)
#' col_order <- unlist(lapply(names(ranked), 
#'                            function(x) grep(paste0("^", x), colnames(logcpm))))
#'                            
#' # use scaled logCPMs                           
#' logcpmZ <- t(scale(t(logcpm[unique(unlist(signatures)),])))
#' pheatmap(mat=logcpmZ[,col_order],
#'          show_rownames=FALSE, cluster_rows=FALSE, cluster_cols=FALSE)
#'          
#' # or each signature individually:
#' pheatmap(mat=t(scale(t(logcpm[signatures$CD4T,]))), show_rownames=FALSE)
#' pheatmap(mat=t(scale(t(logcpm[signatures$CD8T,]))), show_rownames=FALSE)
#' pheatmap(mat=t(scale(t(logcpm[signatures$NveB,]))), show_rownames=FALSE)
#' pheatmap(mat=t(scale(t(logcpm[signatures$NK,]))), show_rownames=FALSE)
#'          
#' # or genes that separate the CD4T and CD8T cells from the rest
#' signatures2 <- CreateGeneSignatures(ranked=ranked, use_groups=c("CD4T", "CD8T"), extended=FALSE)
#' logcpmZ2 <- t(scale(t(logcpm[signatures2,])))
#' pheatmap(mat=logcpmZ2[,col_order],
#'          show_rownames=TRUE, cluster_rows=FALSE, cluster_cols=FALSE)
#'          
#' # find signatures for each group but exclude the CD4T categorically
#' signatures3 <- CreateGeneSignatures(ranked=ranked, exclude_groups=c("CD4T"), extended=FALSE)
#' 
#' # Find all markers for each group against all but one group and output the extended table
#' # to easily see which group the gene does not qualify as a marker for.
#' # Here, the genes in the below tail command are markers for CD4T against NK and NveB but do not separate CD4T from CD8T.
#' signatures4 <- CreateGeneSignatures(ranked=ranked, min.prop=2/3, extended=TRUE)
#' tail(signatures4$CD4T)
#' 
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
#' In version 2.0.0 we introduced the \code{use_groups} argument. This allows to specify a group of celltypes and the function will then infer 
#' markers for this group versus the rest. In the examplle aboves we use it to define CD4/CD8 markers. The advantage here is that these
#' cells have many markers in common so by defining them as a group we only look for "pan"-T cell markers and that retains more genes
#' that when treating both celltypes as independent groups. 
#' 
#' @export
CreateGeneSignatures <- function(ranked, 
                                 use_groups=NULL, exclude_groups=NULL,
                                 delim="_vs_", keep.n=Inf, min.prop=1,
                                 extended=FALSE){
  
  #/ checks
  if(!is.null(use_groups)){
    if(!sum(use_groups %in% names(ranked))==length(use_groups))
      stop("Make sure use_groups are part of names(ranked)")
  }
  
  if(!is.null(exclude_groups)){
    if(!sum(exclude_groups %in% names(ranked))==length(exclude_groups))
      stop("Make sure exclude_groups are part of names(ranked)")
  }
  
  if(!is.null(use_groups) & !is.null(exclude_groups)){
    if(length(intersect(use_groups, exclude_groups))>1)
      stop("There are overlaps between elements in use_groups and exclude_groups!")
  }
  
  if(is.null(names(ranked))) stop("ranked has no names")
  if(min.prop > 1 | min.prop < 0) stop("min.prop must be between > 0 and <= 1")
  
  #/ optionally remove exclude_groups from ranked
  if(!is.null(exclude_groups)){
    
    keep_them <- setdiff(names(ranked), exclude_groups)
    
    ranked <- 
      sapply(keep_them, function(x){
        
        r <- ranked[[x]]
        r[setdiff(keep_them, x)]
        
      }, simplify=FALSE)
    
  }
  
  #/ optionally only look at the group(s) in use_groups rather than all vs all
  if(is.null(use_groups)){
    
    nms <- names(ranked)
    
  } else nms <- use_groups
  
  #/ find the signature genes
  s <- sapply(nms, function(x){
    
    genes.ranged <- ranked[[x]]
    
    remove_these <- if(is.null(use_groups)) x else nms
    
    genes.ranged[remove_these] <- NULL # if the list contains the reference cell type itself
    
    #/ rank matrix, low values means high ranks and viceversa, 1 means not present in that pairwise list: 
    rmat <- rankmatrix(genes.ranged=genes.ranged)
    
    # Get the genes that qualify as markers given the number of groups and the min.prop
    isna <- !is.na(rmat)
    from <- ncol(rmat)
    to   <- floor(min.prop*from)
    to <- if(to==0) 1 else to
    
    med  <- function(y) median(y,na.rm=TRUE)
    
    #/ take the rankmatrix (so each entry is the rank of the gene in the individual ranking),
    #/ and then calculate the median of the ranks -- that is the final ranking metric
    final <- lapply(from:to, function(n){
      
      idx <- base::rowSums(isna)==n
      if(sum(idx)==0) return(NULL)
      
      rmat_x <- rmat[idx,,drop=FALSE]
      
      #/ get the median of the ranks
      mat <- apply(rmat_x, 1, med)
      
      #/ rank by the median
      o <- order(mat)
      
      markers.df <- data.frame(rmat_x[o,,drop=FALSE], check.names=FALSE)
      
      return(markers.df)
      
    })
    
    final <- do.call(rbind, final)
    
    final[!is.na(final)] <- 1
    final[is.na(final)]  <- 0
    final <- head(final, n=keep.n)
    
    return(final)
    
  },simplify = FALSE)
  
  #/ if use_groups is set then only return genes that qualify as markers in all members of 
  #/ use_groups with the same "off-target" pattern so if min.prop demands gene being a marker
  #/ against all but one group then this "one" group must be the same for all elements of use_groups
  if(!is.null(use_groups) & length(use_groups)>1 & min.prop < 1){
    
    #/ strict intersect
    isec <- Reduce(intersect, lapply(as.list(s), function(x) rownames(x)))
    
    if(length(isec)==0) {
      message("No signature genes found!")
      return(NULL)
    }
    
    #/ get the genes
    lp <- sapply(s, function(x) x[isec,], simplify=FALSE)
    cb <- combn(names(s), 2)
    
    lpx_check <- 
      lapply(1:ncol(cb), function(x) {
        
        y <- cb[,x]
        lpx <- lp[[y[1]]]==lp[[y[2]]]
        
        rownames(lpx[rowSums(lpx)==ncol(lpx),])
        
      })
    
    isec_x <- Reduce(intersect, lpx_check)
    
    to_return <- s[[1]][isec_x,]
    
  } else to_return <- s
  
  if(!is.null(use_groups) & min.prop==1 & length(use_groups)>0){
    
    j <- Reduce(intersect, lapply(to_return, rownames))
    to_return <- to_return[[1]][j,]
    
  }
  
  #/ either return the full table or just the names
  if(!extended){
    
    if(class(to_return)[1]=="list"){
      to_return <- sapply(to_return, function(x) rownames(x), simplify=FALSE)  
    } else to_return <- rownames(to_return)
    
    
  }
  
  to_return <- if(length(to_return)==1) to_return[[1]] else to_return
  
  return(to_return)
  
}

#' For the genes being present in genes.ranged per celltype we rank from 1 (best rank) to length 
#' of the current genes.ranged entry, and NA for genes not in that respective entry but in u.
#' Modified from Kolde et al in RobustRankAggreg::rankMatrix()
rankmatrix <- function (genes.ranged) 
{
  u = unique(c(genes.ranged, recursive = TRUE))
  N = length(u)
  
  #/ matrix with only NAs, to be filled with the ranks:
  rmat = matrix(NA, nrow=length(u), ncol=length(genes.ranged), 
                dimnames=list(u, names(genes.ranged)))
  
  #/ fill matrix with ranks per gene and celltype.
  #/ If gene is missing in a celltype (=was not signif vs other celltypes),
  #/ then leave at NA:
  for (i in names(genes.ranged)) {
    
    rmat[genes.ranged[[i]], i] = seq(1, length(genes.ranged[[i]]))
    
  }
  return(rmat)
}

