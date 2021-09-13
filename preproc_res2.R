library(DESeq2)
library(here)

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE177054

dir.create(here("tmp"))

download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE177054&format=file&file=GSE177054%5FGene%5FExpression%2Exlsx",
  destfile = here("tmp/GSE177054.xlslx")
)

counts <- readxl::read_xlsx(here("tmp/GSE177054.xlslx"))

genes <- counts$gene_id

mat <- as.matrix(select(counts, starts_with("count")))
dimnames(mat) <- list(genes, colnames(mat))

meta <- tibble(condition = str_extract(colnames(mat), "DMSO|LEE011|L_D"))

dds <- DESeqDataSetFromMatrix(mat, meta, ~condition)

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "LEE011", "DMSO"))

summary(res)

res <- as_tibble(res, rownames = "gene")

save(res, file = "~/Git/Pathway-Enrichment-Exercises/data/res2.RData")

file.remove(dir(here("tmp"), full.names = TRUE))
file.remove(here("tmp"))
