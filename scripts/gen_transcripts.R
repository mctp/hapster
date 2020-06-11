library(stringr)
library(rtracklayer)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome)
library(foreach)

##
alt.fa <- list.files("refs/hs-hg38/hla/alts", full.names=TRUE, pattern="^[bA-Z].*.fa$")
alt.gff <- list.files("refs/hs-hg38/hla/gff", full.names=TRUE, pattern="^[bA-Z].*.gff$")

genes <- mapply(list, alt.fa, alt.gff, SIMPLIFY=F)
names(genes) <- str_replace(basename(alt.fa), ".fa", "")
foreach(gene=names(genes)) %do% {
    gene.fa <- genes[[gene]][[1]]
    gene.gff <- genes[[gene]][[2]]
    gene.seq <- readDNAStringSet(gene.fa)
    gene.gr <- import(gene.gff)
    gene.exon <- gene.gr[gene.gr$type %in% c("exon", "five_prime_UTR", "three_prime_UTR")]
    gene.exon$type <- "exon"
    gene.txs <- split(gene.exon, seqnames(gene.exon))
    tx.seq <- extractTranscriptSeqs(gene.seq, gene.txs)
    writeXStringSet(tx.seq, sprintf("refs/hs-hg38/hla/rna/%s.fa", gene))
}
