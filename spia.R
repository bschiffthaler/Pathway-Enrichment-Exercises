options(warn = -1)
options(connectionObserver = NULL)

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(tidyverse)
  library(here)
  library(SPIA)
  library(magrittr)
})

load(here("data/res.RData"))

#' Most GSEA packages work best if genes are annotated with entrez identifiers
#' so we first use `AnnotationDbi` to convert our data.
entrezid <- AnnotationDbi::select(org.Hs.eg.db, res$gene, "ENTREZID", "SYMBOL")
res_entrez <- left_join(res, entrezid, by = c("gene" = "SYMBOL"))
#' When converting it is often unavoidable to lost some data due to non 
#' identical annotations
res_entrez %<>% 
  filter(!is.na(ENTREZID) & !duplicated(ENTREZID)) %>%
  filter((!near(baseMean, 0)) & (!is.na(log2FoldChange)))

deg <- res_entrez$log2FoldChange
names(deg) <- res_entrez$ENTREZID

#' Now let's run the SPIA algorithm. Two big caveats: firstly, if we didn't
#' create an up-to-date datastructure for the KEGG annotation, this will use
#' the by now ancient freely available ids, so in order to use this 
#' package effectively a license is more than recommended. To do this, see
#' chapter 3 in the manual. Secondly, it's slow
res_spia <- spia(deg, keys(org.Hs.eg.db))

head(res_spia[, 1:11])
