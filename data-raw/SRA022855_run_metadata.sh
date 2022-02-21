#! /bin/bash

## This code was used to download the metadata for Ravel et. al., 2011 
## (bacteria vaginosis dataset 1). 
## Installation of E-utilities is required. Check:
## https://www.ncbi.nlm.nih.gov/books/NBK179288/
esearch -db sra -query SRA022855 | efetch -format runinfo > SRA022855_run_metadata.csv
