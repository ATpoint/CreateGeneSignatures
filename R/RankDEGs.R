#' Rank genes based on differential expression statistics
#' 
#' Rank genes from lists of pairwise comparison by significance or effect size
#' 
#' @param res a named list of pairwise DE results, see details.
#' @param delim a string that delimits the comparison groups in \code{names(res)}, e.g. celltype1_vs_ccelltype2 would be "_vs_"
#' @param signif.column colname storing significances to use for filtering, e.g. FDR
#' @param signif.threshold keep only genes with \code{signif.column} below that threshold
#' @param effect.column colname storing the effect size, e.g. logFC. Must be a zero-centered effect size so effect size > 0 means higher in one group and < 0 means lower.
#' Don't use something like AUCs from a Wilcox test where > 0.5 means higher and < 0.5 means lower per group.
#' @param effect.threshold keep only genes with effect.column above this threshold, could be a minimum effect size even though it is recommended to explicitely test against
#' the desired minimum effect size rather than postfiltering, see details.
#' @param gene.column colname storing genes or any kind of row identifiers, those will be returned in the output meeting the above criteria
#' @param rnk.method the ranking method for genes meeting criteria, see details.
#' 
#' @details 
#' For an example of how the input should look like see the examples. The pairwise comparisons must be unique, so if something like
#' celltype1_vs_celltype2 is present then do not include celltype2_vs_celltype1 into \code{res} as this is indentical and only the sign of the 
#' effect size changes. The function handles this internally for every celltype.
#' 
#' The ranking has three options:
#' 1) By significance if setting \code{rnk.method="signif"}, using signed \code{-log10(signif.column)} ranking low significances high and vice versa
#' This approach makes sense especially if significances come from a function such as \code{glmTreat} from `edgeR` where a minimum fold change was used
#' as Null hypothesis to eliminate genes with potentially high significances but small effect sizes such as genes with large expression values.
#' See for details the TREAT paper (https://doi.org/10.1093/bioinformatics/btp053) and an answer from the senior author at StackExchange Bioinformatics
#' (https://bioinformatics.stackexchange.com/a/13580/3051). This is only a suggestion, the user is of course free to use any testing machinery.  
#' 2) By effect sizes if setting \code{rnk.method="effect"}, using the \code{-log10(effect.column)} ordered decreasingly. 
#' This approaches probabaly should only be used if the effect sizes were corrected (shrunken) to avoid large effects due to small counts, as the effect size
#' estimates are usually noisy when counts are low. See for more background information e.g. the `DESeq2` vignette towards effect size shrinkage or this paper
#' researching and discussion effect size shrinkage (https://doi.org/10.1093/bioinformatics/bty895).  
#' 3) By a combination of the aforementioned settings if \code{rnk.method="combination"}. In this case the ranking is based on 
#' \code{-log10(effect.column)} * \code{effect.column} ordered decreasingly. This may might sense if one wants to include effect sizes into the ranking
#' but penalize large effect sizes with low significances. 
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
#' For an example with real data see the examples of the \code{CreateGeneSignatures} function of this package.
#' 
#' @author Alexander Toenges
#' 
#' @examples 
#' # first make some example DE results, then rank:
#' res <- sapply(c("gr1_vs_gr2","gr2_vs_gr3","gr1_vs_gr3"), function(x){
#'   data.frame(Gene=paste0("Gene",1:10), logFC=rnorm(10,1,2),FDR=jitter(rep(0.04, 10), 20))
#' },simplify=FALSE)
#' 
#' ranked <- RankDEGs(res = res)
#' 
#' @export
RankDEGs <- function(res, delim="_vs_", 
                     signif.column="FDR", signif.threshold=0.05,
                     effect.column="logFC", effect.threshold=0,
                     gene.column="Gene", rnk.method=c("signif", "effect", "combi")){
  
  ####################################
  # Checks
  ####################################
  
  if(class(res) != "list" | is.null(names(res))){
    stop("res must be a named list", call.=FALSE)
  }
  
  if(!all(grepl(delim, names(res)))) stop("delim was not found in all names of the res")
  
  invisible(match.arg(arg = class(signif.threshold), choices = "numeric"))
  invisible(match.arg(arg = class(effect.threshold), choices = "numeric"))
  
  #/ Fields existing:
  rnk.method <- match.arg(rnk.method)
  
  check.gf <- sum(unlist(lapply(res, function(x) 
    if(!gene.column %in% colnames(x)) return(1) else return(0))))
  check.sf <- sum(unlist(lapply(res, function(x) 
    if(!signif.column %in% colnames(x)) return(1) else return(0))))
  check.ef <- sum(unlist(lapply(res, function(x) 
    if(!effect.column %in% colnames(x)) return(1) else return(0))))
  
  if(check.gf>0) stop("gene.column does not exist in all entries of res.list!")
  if(check.sf>0) stop("signif.column does not exist in all entries of res.list!")
  if(check.ef>0) stop("effect.column does not exist in all entries of res.list!")
  
  ####################################
  # Ranking
  ####################################
  
  #/ make robust against syntactically non-valid names
  unq <- sort(unique(unlist(strsplit(names(res), delim))))
  lookup <- data.frame(original=unq, new=make.names(unq))
  
  if(!all(lookup$original==lookup$new))
    warning(paste("Detected syntactically invalid names.",
                  "It will work anyway, but be sure to double-check your results!",
                  sep="\n"))
  
  unq <- make.names(unq)
  names(res) <- make.names(names(res))
  
  l <- list()
  for(i in unq){
    
    nm <- grep(paste(paste0("^",i,delim),paste0(delim,i,"$"), sep="|"), 
               names(res), value = TRUE)
    
    l[[i]] <- 
      lapply(X = nm, FUN = function(x)
        
      {
        
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
          
          if(rnk.method == "effect"){
            ranked <- tt[order(tt[,effect.column], decreasing=TRUE),][,gene.column]
          }
          if(rnk.method == "signif"){
            ranked <- tt[order(tt[,signif.column], decreasing=FALSE),][,gene.column]
          }
          if(rnk.method == "combi"){
            tt[,"combi"] <- -log10(tt[,signif.column]+.Machine$double.xmin)*tt[,effect.column]
            ranked <- tt[order(tt[,"combi"], decreasing=TRUE),][,gene.column]
          }
          
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
