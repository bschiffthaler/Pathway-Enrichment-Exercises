# Pathway Enrichment Exercises

I have placed a `DESeq` result from [another study](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE177054)
in `data/res2.RData`. Analyze the result with each of ORA, FCS and SPIA methods.

How you do that is relatively freeform, but I suggest you start with the scripts `ora.R`, `fcs.R` and `spia.R`
and modify them to suit your needs. But beware! The gene identifiers used are different so make sure `AnnotationDbi`
gets the right settings to convert what's in `res2.RData` to ENTREZ identifiers.