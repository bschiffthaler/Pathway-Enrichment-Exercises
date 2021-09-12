options(warn = -1)

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(tidyverse)
  library(here)
  library(SPIA)
  library(magrittr)
})

load(here("data/sn_nsm.RData"))

#' Most GSEA packages work best if genes are annotated with entrez identifiers
#' so we first use `AnnotationDbi` to convert our data.
entrezid <- AnnotationDbi::select(org.Hs.eg.db, smokers$ID, "ENTREZID", "ENSEMBL")
smokers_entrez <- left_join(smokers, entrezid, by = c("ID" = "ENSEMBL"))
#' When converting it is often unavoidable to lost some data due to non 
#' identical annotations
smokers_entrez %<>% filter(!is.na(ENTREZID) & !duplicated(ENTREZID))

deg <- smokers_entrez$logFC
names(deg) <- smokers_entrez$ENTREZID

#' Now let's run the SPIA algorithm. Two big caveats: firstly, if we didn't
#' create an up-to-date datastructure for the KEGG annotation, this will use
#' the by now ancient freely available ids, so in order to use this 
#' package effectively a license is more than recommended. To do this, see
#' chapter 3 in the manual. Secondly, it's slow
res_spia <- spia(deg, keys(org.Hs.eg.db))

head(res_spia[, 1:11])
