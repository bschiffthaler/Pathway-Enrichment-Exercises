options(warn = -1)
options(connectionObserver = NULL)

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(KEGGREST)
  library(magrittr)
  library(tidyverse)
  library(magick)
  library(pathview)
  library(fgsea)
  library(here)
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


#' Next, let's extract the feature of interest, in this case the log fold
#' change
scores <- res_entrez$log2FoldChange 
names(scores) <- res_entrez$ENTREZID

#' We can next use `AnnotationDbi` again to create a list of all KEGG 
#' pathways and which genes are represented within them. This will be needed
#' for the next step
pathways <- AnnotationDbi::select(
  org.Hs.eg.db, names(scores),
  "PATH", "ENTREZID"
) %>% as_tibble() %>%
  filter(!is.na(PATH))
pathways <- split(pathways$ENTREZID, pathways$PATH)

#' Finally, we can test using `fgsea`
fgsea_res <- fgsea(pathways = pathways, stats = scores, eps = 0, nproc = 4)
fgsea_res <- arrange(fgsea_res, padj)
head(fgsea_res)

#' `fgsea` has some modest visualization features. Here we plot the leading
#' edge of the most significant pathway (Cell Cycle: has04110)

#+ fig.width=10, fig.height=6
plotEnrichment(pathways[[fgsea_res$pathway[1]]], scores) +
  labs(title = fgsea_res$pathway[1])

#' First we select the top 10 most up and downregulated pathways. We also reverse
#' the order of the downregulated list to start with the least significant,
#' this is mostly for visual aesthetic
top_pathways_up <- fgsea_res[ES > 0][head(order(pval), n=10), ]
top_pathways_down <- fgsea_res[ES < 0][head(order(pval), n=10), ][order(pval, decreasing = TRUE)]
top_pathways <- rbind(top_pathways_up, top_pathways_down)

#' Now we can use KEGGREST to translate the numeric pathway ids into their names.
#' The `keggGet` function actually gets a lot more info, but we'll be content
#' with the name
t_up_names <- keggGet(paste0("hsa", top_pathways_up$pathway))
t_up_names <- sapply(t_up_names, "[[", "NAME") 
t_down_names <- keggGet(paste0("hsa", top_pathways_down$pathway))
t_down_names <- sapply(t_down_names, "[[", "NAME")

pw_names <- sapply(str_split(c(t_up_names, t_down_names), " - "), "[[", 1)

#' We also need to translate the pathWay ids in a number of auxiliary data
#' structures which are required by the plotting function. This is somewhat
#' of a bodge, but the only legal way since KEGG has changed their license
#' terms
pathways_trans <- pathways[top_pathways$pathway]
top_pathways$num_id <- top_pathways$pathway
top_pathways$pathway <- pw_names
names(pathways_trans) <- top_pathways$pathway

#' Finally, we can plot the leading and trailing edges

#+ fig.width=10, fig.height=6
plot.new()
plotGseaTable(
  pathways_trans[top_pathways$pathway], scores, top_pathways, 
  gseaParam=0.5
)

#' If we want detailed info on a pathway, we can make use of the pathview
#' package. Let's zoom in on the p53 (04115) pathway.
pathview(gene.data = scores, pathway.id = "04115",
         low = list("gene" = "skyblue"),
         mid = list(gene = "white"),
         high = list("gene" = "firebrick"))


img <- image_read(here("hsa04115.pathview.png"))
#+ fig.width=10, fig.height=6
img

#' The fgsea package allows us also to quickly get Reactome annotations for
#' our gene list and test those
r_pathways <- reactomePathways(names(scores))
r_fgsea_res <- fgsea(
  pathways = r_pathways, stats = scores, eps = 0, nproc = 4
)
r_fgsea_res <- arrange(r_fgsea_res, padj)
#' We can again plot the leading edge
r_top_pathways_up <- r_fgsea_res[ES > 0][head(order(pval), n=10), ]
r_top_pathways_down <- r_fgsea_res[ES < 0][head(order(pval), n=10), ][
  order(pval, decreasing = TRUE)]
r_top_pathways <- rbind(r_top_pathways_up, r_top_pathways_down)

#+ fig.width=10, fig.height=6
plot.new()
plotGseaTable(
  r_pathways[r_top_pathways$pathway], scores, r_fgsea_res, 
  gseaParam=0.5
)

#' If there is a lot of undesired redundancy in the results, `fgsea` offers
#' functionality to focus on main pathways only
collapsed_pathways <- collapsePathways(r_fgsea_res[order(pval)][padj < 0.01], 
                                       r_pathways, scores)
main_pathways <- r_fgsea_res[pathway %in% collapsed_pathways$mainPathways][
  order(-NES), pathway]

#+ fig.width=10, fig.height=6
plot.new()
plotGseaTable(
  r_pathways[main_pathways], scores, r_fgsea_res, gseaParam=0.5
)
