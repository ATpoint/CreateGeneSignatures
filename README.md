# CreateGeneSignatures

This package implements a simple yet effective rank-based strategy to produce 
celltype-specific gene signatures based on RNA-seq differential expression results.

## Example workflow:
For this example we assume that the package is installed as described below.
We use RNA-seq from [Haemopedia](https://www.haemosphere.org/datasets/show), in this case four lymphoid celltypes from `Haemopedia-Human-RNASeq` dataset,
which are included in the package example data. We first perform differential analysis with `edgeR` and then create a signature for every of the four celltypes,
namely CD4 T cells, CD8 T cells, naive B cells and natural killer (NK) cells.

We first perform differential analysis between all celltypes, then rank the genes for every comparison based on significance, and then derive signatures
based on these lists of ranked genes.

```{r}

# load RNA-seq data for CD4T-, CD8T and naive B cells from Haemopedia:
counts <- readRDS(paste0(
            system.file("extdata",package="CreateGeneSignatures"),
            "/haemopedia_subset.rds"))

# Use edgeR to perform all pairwise comparisons

library(edgeR)
y <- DGEList(counts=counts,group=gsub("\\..", "", colnames(counts)))
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
signatures <- CreateGeneSignatures(ranked=ranked, keep.n=50, min.prop=1, extended=FALSE)
# check number of genes. for CD8T cells we found < 50 genes:
lapply(signatures,length) 


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
         
```

A heatmap of the combined signatures genes:

![heatmap](misc/heatmap.png)


The signatures represent those combination of genes that best separate each of the individual celltypes from all other celltypes. 
In the above example the parameters were very strict, with `min.prop=1` requiring that signature genes ranked highly in every of the initial pairwise conparisons.
This might make sense if celltypes are very different from each otherand one aims to derive a conservative set of marker genes.

In cases where celltypes are less separated (above example represents terminally-differentiated cells), e.g. celltypes that are closely related in
a developmental continuum such as multipotent stem/progenitor populations, such stringent settings might lead to few or no signature genes at all.
In this case lowering the `min.prop` parameter is necessary. As discussed above a signature aims to separate celltypes based on the combination or genes rather
than requiring every single genes to be strictly upregulated in a given celltype versus all other celltypes.

The user is therefore encouraged to try different settings, checking how many genes per celltype can be obtained with the given setting followed by inspection of the separation between celltypes using heatmaps.

## Why ranks?
There are multiple strategies to transform differential expression results into gene signatures. 
Given that differential analysis is the first step in this workflow one might be tempted to simply combine the obtained p-values into a single value, e.g.
using Fisher's method, in order to select genes based on that value. The problem with p-values is that they, beside the effect size, are also a function of statistical power, influenced by e.g. the expression level of a gene and its length. Also the number of samples per group play a notalbe role. 
A gene with identical expression level and effect size measured in a 10 vs 10 sample comparison will result in lower p-values compared to a 3 vs 3 sample comparison. This becomes an even greater problem when dealing with single-cell data if doing single-cell level DE testing, e.g. comparing clusters with 20 and 200 cells versus clusters with 500 and 2000 cells. We therefore chose to only use the DE statistics to rank the differential genes per comparison and then use these ranks for the downstream analysis. This avoids directly comparing p-values between groups with differences in statistical power. The same problem holds true when comparing effect sizes rather than p-values as effect sizes are unreliable in case of low replicate numbers and/or the presence of low counts.

## Installation

```r

# The package is fully base R with no dependencies.

install.packages("remotes")
remotes::install_github("ATpoint/CreateGeneSignatures")

```
