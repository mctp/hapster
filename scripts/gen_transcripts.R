library(stringr)
library(rtracklayer)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome)
library(foreach)

alt.fa <- snakemake@input[["fastas"]]
alt.gff <- snakemake@input[["gffs"]]
gene <- str_replace(basename(alt.fa), ".fa", "")

seq <- readDNAStringSet(alt.fa)
gr <- import(alt.gff)
exon <- gr[gr$type %in% c("exon", "five_prime_UTR", "three_prime_UTR")]
exon$type <- "exon"
txs <- split(exon, seqnames(exon))
txs.seq <- extractTranscriptSeqs(seq, txs)
writeXStringSet(txs.seq, sprintf("refs/hs-hg38/hla/rna/%s.fa", gene))
