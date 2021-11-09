library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

#Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
total_rna_filename <- args[1]
hla_rna_filename <- args[2]
output_filename <- args[3]

total_rna <- read_tsv(total_rna_filename)
hla_rna <- read_tsv(hla_rna_filename)

recalc_rna <- bind_rows(total_rna, hla_rna) %>%
  mutate(rate = est_counts/eff_length,
         tpm_recalc = (est_counts/eff_length)/sum(rate)*1000000) %>%
  dplyr::select(target_id, length, eff_length, est_counts, tpm_recalc) %>%
  dplyr::rename("tpm" = tpm_recalc) %>%
  filter(grepl("\\*", target_id))

write_csv(recalc_rna, output_filename)