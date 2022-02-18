# Analyses using MicrobiomeBenchmarkData 

This package contains analyses performed on datasets included or under
consideration to be included in the 
[MicrobiomeBenchmarkData](https://waldronlab.io/MicrobiomeBenchmarkData)
package.

```{r}
## Package installation
if (!"BiocManager" %in% install.packages()[,"package"])
  install.packages("BiocManager")
BiocManager::install("waldronlab/MicrobiomeBenchmarkDataAnalyses")
```