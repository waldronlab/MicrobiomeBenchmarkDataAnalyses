library(ZINQ)
library(MicrobiomeBenchmarkData)
library(magrittr)


data(Sample_Data)

tax_tab <- Sample_Data['CSS_taxon2']
metadata <- Sample_Data['X']

check_output <- ZINQ_check(tax_tab = tax_table, metadata = metadata, C = 'X')
