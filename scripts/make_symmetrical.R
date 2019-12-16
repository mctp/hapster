library(readr)
library(dplyr)

matrix_raw_filename <- snakemake@input[["matrix_raw"]]
matrix_filename <- snakemake@output[["matrix"]]

make_sym <- function(nm, gene) {
  likelihoods <- read_csv(matrix_raw_filename)
  likelihoods <- likelihoods[, 2:ncol(likelihoods)]
  
  likelihoods_values <- as.matrix(likelihoods)
  likelihoods_values_t <- t(likelihoods_values)
  
  likelihoods_values_mean <- (likelihoods_values + likelihoods_values_t) / 2
  
  likelihoods_mean <- data.frame(x = names(likelihoods), stringsAsFactors = FALSE) %>%
    cbind(likelihoods_values_mean) %>%
    rbind(c("", names(likelihoods)), ., stringsAsFactors = FALSE)
  
  write_csv(likelihoods_mean, matrix_filename,
            col_names = FALSE, quote_escape = FALSE)
}

make_sym(matrix_raw_filename, matrix_filename)
