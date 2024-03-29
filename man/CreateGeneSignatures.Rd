% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreateGeneSignatures.R
\name{CreateGeneSignatures}
\alias{CreateGeneSignatures}
\title{Create gene signatures in a rank-based style}
\usage{
CreateGeneSignatures(
  ranked,
  use_groups = NULL,
  exclude_groups = NULL,
  delim = "_vs_",
  keep.n = Inf,
  min.prop = 1,
  extended = FALSE
)
}
\arguments{
\item{ranked}{nested lists of ranked genes for each celltype, see details and the RankDEGs function}

\item{use_groups}{trigger group-aware mode. A vector with celltype names to be used to define markers for as a group,
so this group versus all others. See details.}

\item{exclude_groups}{vector, remove these groups from the signature creation}

\item{delim}{a string indicating the delimiter of names(ranked), e.g. celltype1_vs_celltype2 would be "_vs_"}

\item{keep.n}{number of genes to keep per signature, by default all candidates are returned for manual posthoc filtering}

\item{min.prop}{minimum proportion of comparisons per celltype that a gene must be included in. See details.}

\item{extended}{logical, whether to output the signature genes as a data.frame that indicates whether the gene
qualified as a marker against the other groups (1) or not (0)}
}
\description{
Convert lists of ranked genes into celltype-specific signatures
}
\details{
The function takes as input a nested list that, for every celltype, contains the ranked genes from pairwise comparisons against every other
celltype in this analysis, see also the \code{RankDEGs} function of this package. For every celltype first a rank matrix is created that
stores the ranked genes from every pairwise comparison, assigning low scores to genes that ranked high in the input ranking and vice versa.
Genes to form the signatures are then selected by choosing those genes with minimal ranks across all comparisons.
In an iterative fashion genes that ranked high all comparisons are chosen first (and then ordered by median score), then genes ranked high in all but one, 
in all but two, all but three etc comparisons. The user can choose via \code{min.prop} the fraction of comparisons that a gene be highly ranked in.
For example a setup with four celltypes and \code{min.prop=1} would mean that only genes are included that ranked high in all comparisons.
If \code{min.prop=2/3} is chosen then signature genes must rank high in at least two out of three comparisons. The higher this value between (between zero and one)
the more stringent the filtering. If celltypes are decently separated it might be trivial to find large numbers of signature genes even with min.prop of 1.
In cases where cells are very similar to each other, e.g. cells being close in a developmental continuum, one will need to lower this threshold in order to get some
(or any) signature genes. 

The \code{use_groups} argument allows to specify a group of celltypes and the function will then infer 
markers for this combined group versus the rest. In the example above we use it to define combined CD4/CD8 markers. The advantage here is that these
cells have many markers in common so by defining them as a group we only look for "pan"-T cell markers and that retains more genes
that when treating both celltypes as independent groups.
}
\examples{
# load RNA-seq data for CD4T-, CD8T and naive B cells from Haemopedia:
counts <- readRDS(paste0(
            system.file("extdata",package="CreateGeneSignatures"),
            "/haemopedia_subset.rds"))

# Use edgeR to perform all pairwise comparisons

library(edgeR)
y <- DGEList(counts=counts,group=gsub("\\\\..", "", colnames(counts)))
design <- model.matrix(~0+group,y$samples)
colnames(design) <- gsub("group", "", colnames(design))
y <- y[filterByExpr(y),,keep.lib.size=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)

# all unique pairwise contrasts:
contrasts <- makeContrasts(CD4T_vs_CD8T  = CD4T-CD8T,
                           CD4T_vs_NK    = CD4T-NK,
                           CD4T_vs_NveB  = CD4T-NveB,
                           CD8T_vs_NK    = CD8T-NK,
                           CD8T_vs_NveB  = CD8T-NveB,
                           NK_vs_NveB    = NK-NveB,
                           levels = design)
                           
# test using glmTreat against a minumum fold change of ~1.5                          
res <- sapply(colnames(contrasts), function(con){
  tt<-topTags(glmTreat(fit,contrast=contrasts[,con],lfc=log2(1.5)),n=Inf)$table
return(data.frame(Gene=rownames(tt), tt))
}, simplify = FALSE)

# Rank the DEGs:
ranked <- RankDEGs(res)

# Create signatures, keeping top 50 signature genes that separate the respective celltype
# from all other celltypes:
signatures <- CreateGeneSignatures(ranked=ranked, keep.n=50)
# check number of genes. for CD8T cells we found < 50 genes:
lengths(signatures)

# Inspect signatures using heatmaps plotting the scaled logcpms of the signature genes
library(pheatmap)
logcpm <- log2(edgeR::cpm(y,log=FALSE)+1)

# plot a heatmap in the order of names(signatures)
col_order <- unlist(lapply(names(ranked), 
                           function(x) grep(paste0("^", x), colnames(logcpm))))
                           
# use scaled logCPMs                           
logcpmZ <- t(scale(t(logcpm[unique(unlist(signatures)),])))
pheatmap(mat=logcpmZ[,col_order],
         show_rownames=FALSE, cluster_rows=FALSE, cluster_cols=FALSE)
         
# or each signature individually:
pheatmap(mat=t(scale(t(logcpm[signatures$CD4T,]))), show_rownames=FALSE)
pheatmap(mat=t(scale(t(logcpm[signatures$CD8T,]))), show_rownames=FALSE)
pheatmap(mat=t(scale(t(logcpm[signatures$NveB,]))), show_rownames=FALSE)
pheatmap(mat=t(scale(t(logcpm[signatures$NK,]))), show_rownames=FALSE)
         
# or genes that separate the CD4T and CD8T cells from the rest
signatures2 <- CreateGeneSignatures(ranked=ranked, use_groups=c("CD4T", "CD8T"), extended=FALSE)
logcpmZ2 <- t(scale(t(logcpm[signatures2,])))
pheatmap(mat=logcpmZ2[,col_order],
         show_rownames=TRUE, cluster_rows=FALSE, cluster_cols=FALSE)
         
# find signatures for each group but exclude the CD4T categorically
signatures3 <- CreateGeneSignatures(ranked=ranked, exclude_groups=c("CD4T"), extended=FALSE)

# Find all markers for each group against all but one group and output the extended table
# to easily see which group the gene does not qualify as a marker for.
# Here, the genes in the below tail command are markers for CD4T against NK and NveB but do not separate CD4T from CD8T.
signatures4 <- CreateGeneSignatures(ranked=ranked, min.prop=2/3, extended=TRUE)
tail(signatures4$CD4T)

}
\author{
Alexander Toenges
}
