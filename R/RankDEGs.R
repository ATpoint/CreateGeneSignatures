#' Rank genes based on differential expression statistics
#' 
#' Rank genes from lists of pairwise comparison by significance or effect size
#' 
#' @param res a named list of pairwise DE results, see details.
#' @param delim a string that delimits the comparison groups in \code{names(res)}, e.g. celltype1_vs_celltype2 would be "_vs_"
#' @param signif.column colname storing significances to use for filtering, e.g. FDR
#' @param signif.threshold keep only genes with \code{signif.column} below that threshold
#' @param effect.column colname storing the effect size, e.g. logFC. Must be a zero-centered effect size so effect size > 0 means higher in one group and < 0 means lower.
#' Don't use something like AUCs from a Wilcox test where > 0.5 means higher and < 0.5 means lower per group.
#' @param effect.threshold keep only genes with effect.column above this threshold, could be a minimum effect size even though it is recommended to explicitely test against
#' the desired minimum effect size rather than postfiltering, see details.
#' @param gene.column colname storing genes or any kind of row identifiers, those will be returned in the output meeting the above criteria
#' @param rnk.column use this column for the ranking
#' @param rnk.method either 'decreasing' or 'increasing' ranking based on \code{rnk.column}
#' 
#' @details 
#' For an example of how the input should look like see the examples. The pairwise comparisons must be unique, so if something like
#' celltype1_vs_celltype2 is present then do not include celltype2_vs_celltype1 into \code{res} as this is identical and only the sign of the 
#' effect size (e.g. the logFC) changes. The function handles this internally for every celltype.
#' 
#' The \code{signif.column} and \code{effect.column} are used first to filter the data, e.g. for FDR and logFC, and then the ranking is done based
#' on the \code{rnk.column} but aware of the direction of change based on \code{effect.column}, and in its current state only genes
#' with a positive effect size are taken into account. The \code{rnk.column} could e.g. e the nominal PValue or t-stat column which both
#' have the advantage over FDR that they usually have no ties.
#'
#' The output will be a nested list with the ranked genes for every celltype compared to every other celltype based in the entries of \code{res},
#' see the examples. Something like:
#' $celltype1
#' ..$celltype2
#' ..$celltype3
#' $celltype2
#' ..$celltype1
#' ..$celltype3
#' $celltype3
#' ..celltype1
#' ..celltype2
#'
#' For an example with real data see the examples of the \code{CreateGeneSignatures} function of this package.#' 
#' @author Alexander Toenges
#' 
#' @examples 
#' # first make some dummy DE results, then rank:
#' set.seed(1)
#' res <- sapply(c("gr1_vs_gr2","gr2_vs_gr3","gr1_vs_gr3"), function(x){
#'   data.frame(Gene=paste0("Gene",1:10), 
#'              logFC=rnorm(10,1,2),
#'              PValue=jitter(rep(0.001, 10), 20),
#'              FDR=jitter(rep(0.04, 10), 20))
#' },simplify=FALSE)
#' 
#' # this is how the results tables look:
#' res$gr1_vs_gr2
#' 
#' ranked <- RankDEGs(res=res, rnk.column="PValue", rnk.method="increasing")
#' 
#' @export
RankDEGs <- function(res, delim="_vs_", 
                     signif.column="FDR", signif.threshold=0.05,
                     effect.column="logFC", effect.threshold=0,
                     gene.column="Gene", rnk.column="PValue", rnk.method="increasing"){
  
  #---------------------------
  # Checks
  #---------------------------
  
  if(!class(res) %in% c("list", "SimpleList") | is.null(names(res))){
    stop("res must be a named list", call.=FALSE)
  }
  
  if(!all(grepl(delim, names(res)))) stop("delim was not found in all names of the res")
  
  invisible(match.arg(arg=class(signif.threshold), choices="numeric"))
  invisible(match.arg(arg=class(effect.threshold), choices="numeric"))
  
  check.gf <- sum(unlist(lapply(res, function(x) 
    if(!gene.column %in% colnames(x)) return(1) else return(0))))
  check.sf <- sum(unlist(lapply(res, function(x) 
    if(!signif.column %in% colnames(x)) return(1) else return(0))))
  check.ef <- sum(unlist(lapply(res, function(x) 
    if(!effect.column %in% colnames(x)) return(1) else return(0))))
  check.rk <- sum(unlist(lapply(res, function(x) 
    if(!rnk.column %in% colnames(x)) return(1) else return(0))))
  
  if(check.gf>0) stop("gene.column does not exist in all entries of res.list!")
  if(check.sf>0) stop("signif.column does not exist in all entries of res.list!")
  if(check.ef>0) stop("effect.column does not exist in all entries of res.list!")
  if(check.rk>0) stop("rnk.column does not exist in all entries of res.list!")
  
  if(!rnk.method %in% c("decreasing", "increasing"))
    stop("rnk.method must be one of <increasing, decreasing>")
  
  #---------------------------
  # Ranking
  #---------------------------
  
  #/ make robust against syntactically non-valid names
  unq <- sort(unique(unlist(strsplit(names(res), delim))))
  lookup <- data.frame(original=unq, new=make.names(unq))
  
  if(!all(lookup$original==lookup$new))
    warning(paste("Detected syntactically invalid names.",
                  "It will work anyway, but be sure to double-check your results!",
                  "The invalid names are:",
                  paste0(lookup$original[which(!lookup$original==lookup$new)], collapse="\n"),
                  sep="\n"))
  
  unq <- make.names(unq)
  names(res) <- make.names(names(res))
  
  l <- list()
  for(i in unq){
    
    nm <- grep(paste(paste0("^",i,delim),paste0(delim,i,"$"), sep="|"), 
               names(res), value=TRUE)
    
    l[[i]] <- 
      lapply(X=nm, FUN=function(x){
        
        s  <- strsplit(x,delim)[[1]]
        
        tt <- res[[x]]
        
        if(i==s[2]) tt[,effect.column] <- tt[,effect.column]*-1
        
        if(!is.null(signif.column)) {
          sig <- (tt[,signif.column] < signif.threshold)
        } else sig <- rep(TRUE, nrow(tt))
        
        if(!is.null(effect.column)) {
          eff <- (tt[,effect.column] > effect.threshold)
        } else eff <- rep(TRUE, nrow(tt))
        
        tt<-tt[rowSums(data.frame(sig,eff))==2,]
        
        if(nrow(tt)>0){
          
          decreasing <- if(rnk.method=="decreasing") TRUE else FALSE
          
          ix <- sort(tt[,rnk.column], index.return=TRUE)$ix
          
          ranked <- tt[ix,gene.column]
          
        } else ranked <- NULL
        
        return(ranked)
        
      })
    
    names(l[[i]]) <- setdiff(unlist(strsplit(nm,delim)), i)
    
  }
  
  iszero <- sum(!lengths(l) > 0)
  if(iszero>0) message("There are ", iszero, " comparisons with zero genes")
  
  #/ map names back to original
  names(l) <- lookup$original[match(names(l), lookup$new)]
  sapply(l, function(x)  {
    
    names(x) <- lookup$original[match(names(x), lookup$new)]
    x
    
  }, simplify=FALSE)
  
}