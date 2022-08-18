suppressPackageStartupMessages({
  library(eisaR)
  library(Biostrings)
  library(BSgenome)
  library(stringr)
  library(GenomicFeatures)
})

gtf = commandArgs(trailingOnly=TRUE)[1] 

gtf_res = gsub(pattern = ".gtf", "", gtf)
gtf_res = paste0(gtf_res, "_filteredsmRNA.gtf")

hoge = import(con = gtf, format = "gtf") # Import the original exonic genome annotation file
hoge = as.data.frame(hoge)
#remove small RNA annotations
hoge = hoge[!hoge$gene_biotype %in% c("miRNA","snRNA", "snoRNA","scaRNA", "scRNA",  "sRNA", "ribozyme", "vault_RNA", "Mt_tRNA", "Mt_rRNA"), ]

#modified for compatibility to RSeQC
hoge1 = hoge$transcript_id
hoge1[is.na(hoge1)] = ""
hoge$transcript_id = hoge1
hoge = as.data.frame(hoge)
unique(hoge$source)

hoge = makeGRangesFromDataFrame(hoge, keep.extra.columns=TRUE)
rtracklayer::export(hoge, gtf_res, format = "gtf")