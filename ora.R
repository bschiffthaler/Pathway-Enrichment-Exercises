options(warn = -1)

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(tidyverse)
  library(here)
  library(clusterProfiler)
  library(DOSE)
  library(enrichplot)
  library(ReactomePA)
})

load(here("data/sn_nsm.RData"))

#' Most GSEA packages work best if genes are annotated with entrez identifiers
#' so we first use `AnnotationDbi` to convert our data.
entrezid <- AnnotationDbi::select(org.Hs.eg.db, smokers$ID, "ENTREZID", "ENSEMBL")
smokers_entrez <- left_join(smokers, entrezid, by = c("ID" = "ENSEMBL"))
#' When converting it is often unavoidable to lost some data due to non 
#' identical annotations
smokers_entrez %<>% filter(!is.na(ENTREZID) & !duplicated(ENTREZID))

#' We can check which pathways are enriched in our data using the `enrichKEGG`
#' function
kegg_res <- enrichKEGG(smokers_entrez$ENTREZID)
head(kegg_res)

#' If we want detailed info on a pathway, we can make use of the pathview
#' package. Let's zoom in on the strange Malaria (hsa05144) pathway.
gene_data <- smokers_entrez$logFC
names(gene_data) <- smokers_entrez$ENTREZID
pw <- pathview(gene.data = gene_data,
               pathway.id = "hsa05144",
               low = list("gene" = "skyblue"),
               mid = list(gene = "white"),
               high = list("gene" = "firebrick"))
img <- image_read(here("hsa05144.pathview.png"))
#+ fig.width=10, fig.height=6
img

#+ fig.width=10, fig.height=6
enrichplot::dotplot(kegg_res)

#' In some occasions using KEGG modules (which are functional modules) can
#' make results more interpretable
mkegg_res <- enrichMKEGG(smokers_entrez$ENTREZID)
head(mkegg_res)

#' For disease focused research, one sometime wants to explore disease
#' ontology
do_res <- enrichDO(smokers_entrez$ENTREZID)
head(do_res)
enrichplot::dotplot(do_res)

#' Finally, we can also check reactome for overrepresented pathways
pw_res <- enrichPathway(smokers_entrez$ENTREZID)
head(pw_res)
viewPathway("Class B/2 (Secretin family receptors)", 
            readable = TRUE, 
            foldChange = gene_data)
