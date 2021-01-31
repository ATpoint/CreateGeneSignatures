# GeneSignatures

This package implements a simple yet effective rank-based strategy to produce 
celltype-specific gene signatures based on RNA-seq differential expression results.

## Example workflow:
For this example we assume that the package is installed as described below.
We use RNA-seq from (Haemopedia)[https://www.haemosphere.org/datasets/show], in this case four lymphoid celltypes from `Haemopedia-Human-RNASeq` dataset,
which are included in the package example data. We first perform differential analysis with `edgeR` and then create a signature for every of the four celltypes,
namely CD4 T cells, CD8 T cells, naive B cells and natural killer (NK) cells.

We first perform differential analysis between all celltypes, then rank the genes for every comparison based on significance, and then derive signatures
based on these lists of ranked genes.

```{r}

# Low the raw counts:
counts <- readRDS(paste0(system.file("extdata",package="CreateGeneSignatures"), "/haemopedia_subset.rds"))

# Use edgeR to perform all pairwise comparisons
library(edgeR)
y <- DGEList(counts=counts,group=gsub("\\..", "", colnames(counts)))
design <- model.matrix(~0+group,y$samples)
colnames(design) <- gsub("group", "", colnames(design))
y <- y[filterByExpr(y),,keep.lib.size=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)

# define all unique pairwise contrasts:
contrasts <- makeContrasts(CD4T_vs_CD8T  = CD4T-CD8T,
                           CD4T_vs_NK    = CD4T-NK,
                           CD4T_vs_NveB  = CD4T-NveB,
                           CD8T_vs_NK    = CD8T-NK,
                           CD8T_vs_NveB  = CD8T-NveB,
                           NK_vs_NveB    = NK-NveB,
                           levels = design)
                           
# test every contrast against a fold change of 1.5 with glmTreat:
res <- sapply(colnames(contrasts), function(con){

  tt <- topTags(glmTreat(fit, contrast=contrasts[,con], lfc=log2(1.5)), n=Inf)$table
  return(data.frame(Gene=rownames(tt), tt))

}, simplify = FALSE)

# Rank the DEGs:
ranked <- RankDEGs(res)

# Create signatures, keeping the top-50 genes per celltype:
signatures <- CreateGeneSignatures(ranked=ranked, keep.n=50, min.prop=1)

# As can be seen the CD8T cells only have 33 signature genes with the chosen parameters:
lapply(signatures,length) 

# Inspect the signature genes using heatmaps, using the scaled CPMs on the log scale:
library(pheatmap)
logcpm <- log2(edgeR::cpm(y,log=FALSE)+1)
pheatmap(mat=t(scale(t(logcpm[unique(unlist(signatures)),]))), show_rownames=FALSE, cluster_rows = FALSE)

# Alternatively check every signature separately:
pheatmap(mat=t(scale(t(logcpm[signatures$CD4T,]))), show_rownames=FALSE)
pheatmap(mat=t(scale(t(logcpm[signatures$CD8T,]))), show_rownames=FALSE)
pheatmap(mat=t(scale(t(logcpm[signatures$NveB,]))), show_rownames=FALSE)
pheatmap(mat=t(scale(t(logcpm[signatures$NK,]))), show_rownames=FALSE)

```

A heatmap of the combined signatures genes:

![Heatmap](https://i.ibb.co/9VxkTwh/heatmap.jpg)

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
using Fisher's method, in order to select genes based on that value. The problem with p-values is that they, beside the effect size, also a function of statistical power.
A gene with identical expression level and effect size measured in a 10 vs 10 comparison will result in lower p-values compared to a 3 vs 3 comparison.
This becomes an even greater problem when dealing with single-cell data when doing single-cell level DE testing, e.g. comparing clusters with 20 and 200 cells versus
clusters with 500 and 2000 cells. We therefore chose to only use the DE statistics to rank the differential genes per comparison and then use these ranks for the downstream analysis. This avoids directly comparing p-values between groups with differences in statistical power. The same problem holds true when comparing effect sizes
rather than p-values as effect sizes are unreliable in case of low replicate numbers and/or the presence of low counts.

## Installation

```r

# The package is fully base R with no dependencies.

install.packages("remotes")
remotes::install_github("ATpoint/CreateGeneSignatures")

```
