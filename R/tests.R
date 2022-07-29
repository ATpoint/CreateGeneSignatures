install.packages("remotes", repos="https://cloud.r-project.org", verbose=FALSE)
remotes::install_github("atpoint/CreateGeneSignatures")
library(CreateGeneSignatures)

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
