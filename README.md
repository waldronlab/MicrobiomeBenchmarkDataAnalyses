# Analyses using MicrobiomeBenchmarkData 

This package contains analyses performed on datasets included or under
consideration to be included in the 
[MicrobiomeBenchmarkData](https://waldronlab.io/MicrobiomeBenchmarkData)
package, which provides benchmark datasets for differential abundance
methods in microbiome research.

Package website: https://waldronlab.io/MicrobiomeBenchmarkDataAnalyses

Installation

```{r}
## Package installation
if (!"BiocManager" %in% install.packages()[,"Package"])
  install.packages("BiocManager")
BiocManager::install("waldronlab/MicrobiomeBenchmarkDataAnalyses")
```

This package uses the renv package, which could be used to install the
same package versions used here. If that's the case is best to clone
this repository and follow the instructions provided at the 
[renv website](https://rstudio.github.io/renv/articles/renv.html).
