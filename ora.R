options(warn = -1)

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(tidyverse)
  library(here)
  library(clusterProfiler)
  library(DOSE)
  library(enrichplot)
  library(ReactomePA)
  library(pathview)
  library(magick)
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

#' We can check which pathways are enriched in our data using the `enrichKEGG`
#' function
kegg_res <- enrichKEGG(
  filter(res_entrez, !is.na(padj) & padj < 0.05)$ENTREZID,
  universe = res_entrez$ENTREZID
)
head(kegg_res)

#' If we want detailed info on a pathway, we can make use of the pathview
#' package. Let's zoom in on the Spliceosome (hsa03040) pathway.
gene_data <- res_entrez$log2FoldChange
names(gene_data) <- res_entrez$ENTREZID
pw <- pathview(gene.data = gene_data,
               pathway.id = "hsa03040",
               low = list("gene" = "skyblue"),
               mid = list(gene = "white"),
               high = list("gene" = "firebrick"))
img <- image_read(here("hsa03040.pathview.png"))
#+ fig.width=10, fig.height=6
img

#+ fig.width=10, fig.height=6
enrichplot::dotplot(kegg_res)

#' In some occasions using KEGG modules (which are functional modules) can
#' make results more interpretable
mkegg_res <- enrichMKEGG(
  filter(res_entrez, !is.na(padj) & padj < 0.05)$ENTREZID,
  universe = res_entrez$ENTREZID
)
head(mkegg_res)

#' For disease focused research, one sometime wants to explore disease
#' ontology
do_res <- enrichDO(
  filter(res_entrez, !is.na(padj) & padj < 0.05)$ENTREZID,
  universe = res_entrez$ENTREZID
)
head(do_res)

#+ fig.width=10, fig.height=6
enrichplot::dotplot(do_res)

#' Finally, we can also check reactome for overrepresented pathways
pw_res <- enrichPathway(
  filter(res_entrez, !is.na(padj) & padj < 0.05)$ENTREZID,
  universe = res_entrez$ENTREZID
)
head(pw_res)

#+ fig.width=10, fig.height=6
viewPathway("TP53 Regulates Transcription of Cell Cycle Genes", 
            readable = TRUE, 
            foldChange = gene_data)
