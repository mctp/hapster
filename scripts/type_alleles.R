library(readr)
library(tibble)
library(dplyr)

likelihoods <- snakemake@input[["likelihoods"]]
counts <- snakemake@input[["counts"]]
haplotype_filename <- snakemake@output[["haplotype"]]
cutoff <- snakemake@params[['cutoff']]

type_alleles <- function(likelihoods, counts, haplotype_filename, cutoff) {
  likelihoods <-  read_csv(likelihoods)
  likelihoods <- likelihoods[, 2:ncol(likelihoods)]
  counts <-  read_csv(counts)
  
  likelihoods_pared <- as.matrix(likelihoods)
  counts_pared <- counts$count
  
  #Matrix starts too singular, so remove worst offenders (most highly correlated alleles) until it is solvable
  boolFalse<-F
  while(boolFalse==F)
  {
    tryCatch({
      cor_pared <- cor(likelihoods_pared)
      cor_pared[lower.tri(cor_pared, diag = TRUE)] <- 0
      max_cor_pared <- apply(cor_pared, 2, max)
      
      sorted_max <- sort(max_cor_pared, decreasing = TRUE)[1]
      max_column <- c(1:ncol(likelihoods_pared))[colnames(likelihoods_pared) %in% names(sorted_max)]
      max_row <- c(1:nrow(cor_pared))[cor_pared[, max_column] == sorted_max][1]
      
      print(colnames(likelihoods_pared)[max_column])
      likelihoods_pared <- likelihoods_pared[-max_column, -max_column]
      counts_pared <- counts_pared[-max_column]
      
      solve(likelihoods_pared, counts_pared)
      boolFalse<-T
    },error=function(e){
    },finally={})
  }
  
  #Remove high correlation alleles iteratively until only none are left
  denoised_reads <- solve(likelihoods_pared, counts_pared)
  plot(denoised_reads)
  #Get upper triangular matrix of all pairwise correlations
  cor_pared <- cor(likelihoods_pared)
  cor_pared[lower.tri(cor_pared, diag = TRUE)] <- 0
  
  #Create a list of indexes for all pairs with correlation above the cutoff,
  #and find the magnitude of the difference between the two from the denoised reads vector
  indexes <- which(cor_pared > cutoff, arr.ind=TRUE)
  
  while (nrow(indexes) > 1) {
    indexes <- tibble(allele = rownames(indexes),
                      row = indexes[, 1],
                      col = indexes[, 2]) %>%
      mutate(magnitude=abs(denoised_reads[row]-denoised_reads[col]))
    
    #Take the pair with the highest magnitude, and select the member of the pair with the lowest value
    max_mag <- filter(indexes, magnitude==max(magnitude))
    most_negative <- names(denoised_reads[denoised_reads == min(denoised_reads[max_mag$row], denoised_reads[max_mag$col])])[1]
    print(most_negative)
    
    #Remove most negative value
    neg_index <- c(1:length(colnames(likelihoods_pared)))[colnames(likelihoods_pared) %in% most_negative]
    likelihoods_pared <- likelihoods_pared[-neg_index, -neg_index]
    counts_pared <- counts_pared[-neg_index]
    
    #There is a chance the matrix could end up exactly singular, so we skip over that if it errors out
    #If this happens we either use continue on, or use the last valid set of denoised reads
    #if there are no future valid sets
    tryCatch({
      denoised_reads <- solve(likelihoods_pared, counts_pared)
    },error=function(e){
    })
    
    plot(denoised_reads)
    
    #Get upper triangular matrix of all pairwise correlations
    cor_pared <- cor(likelihoods_pared)
    cor_pared[lower.tri(cor_pared, diag = TRUE)] <- 0
    
    #Create a list of indexes for all pairs with correlation above .99,
    #and find the magnitude of the difference between the two from the denoised reads vector
    indexes <- which(cor_pared > cutoff, arr.ind=TRUE)
  }
  
  #Because we perform realignment, we do not make a homozygous/heterozygous call.
  #In the homozygous case, reads will simply realign to the proper allele and mutations can
  #be called.
  denoised_reads <- sort(denoised_reads, decreasing = TRUE)[1:2]
  haplotype <- names(denoised_reads)

  readr::write_csv(tibble(hap = haplotype) %>% arrange(hap), haplotype_filename, col_names = FALSE)
}

type_alleles(likelihoods, counts, haplotype_filename, cutoff)