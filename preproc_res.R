library(DESeq2)
library(here)

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159151

dir.create(here("tmp"))

download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE159151&format=file",
  destfile = here("tmp/GSE159151_RAW.tar")
)

files <- untar(here("tmp/GSE159151_RAW.tar"), exdir = here("tmp"))

files <- dir(here("tmp"), pattern = "*.txt.gz$", full.names = TRUE)

counts <- lapply(files, read_tsv)

genes <- counts[[1]]$`Gene symbol`

counts <- lapply(counts, function(f) {f[,2]}) %>% bind_cols()

mat <- as.matrix(counts)
dimnames(mat) <- list(genes, colnames(counts))

meta <- tibble(condition = sapply(str_split(colnames(mat), "-"), "[[", 1))


dds <- DESeqDataSetFromMatrix(mat, meta, ~condition)

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "plus", "minus"))

res <- as_tibble(res, rownames = "gene")

save(res, file = "~/Git/Pathway-Enrichment-Exercises/data/res.RData")

file.remove(dir(here("tmp"), full.names = TRUE))
file.remove(here("tmp"))