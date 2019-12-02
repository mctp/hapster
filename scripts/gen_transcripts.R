library(rtracklayer)
library(Biostrings)
library(GenomicFeatures)

##
alt.fa <- list.files("refs/hs-hg38/hla/alts", "^[A-Z]", full.names=TRUE)
alt.gff <- list.files("refs/hs-hg38/hla/gff", "^[A-Z]", full.names=TRUE)
alt.seq <- unlist(DNAStringSetList(lapply(alt.fa, readDNAStringSet)))
alt.gr <- unlist(GRangesList(lapply(alt.gff, import)))
alt.exon <- alt.gr[alt.gr$type %in% c("exon", "five_prime_UTR", "three_prime_UTR")]
alt.exon$type <- "exon"
alt.txs <- split(alt.exon, seqnames(alt.exon))
extr <- extractTranscriptSeqs(alt.seq, alt.txs)
writeXStringSet(extr, "refs/hs-hg38/hla/transcripts.fa")
