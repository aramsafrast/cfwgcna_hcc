---
title: "cfRNA & exoRNA module preservation analysis in scRNA data"
subtitle: "Data of Zhu et al. https://doi.org/10.7150%2Fthno.48206, Li et al. https://doi.org/10.1093/nar/gkx891 & MacParland et al. https://doi.org/10.1038/s41467-018-06318-7"
author: "Aram Safrastyan"
date: "23 September, 2022"
output:
  html_document: 
    keep_md: yes
    toc: TRUE
    code_folding: hide
    number_sections: TRUE
---
            


<style type="text/css">
.main-container {
  max-width: 1500px;
  margin-left: auto;
  margin-right: auto;
}
</style>



# load libraries
***


```r
library(Seurat)
library(forcats)
library(dplyr)
library(WGCNA)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(tibble)
library(ggrepel)
library(tidyr)
library(stringr)
library(cowplot)
```

# Module presevation calculation (cfRNA modules => cell specific scRNA data)
## warning - resource intensive
***


```r
enableWGCNAThreads(nThreads=snakemake@threads) 
```

```
## Allowing parallel execution with up to 15 working processes.
```

```r
dir.create("./scrna/results/cf_mp/", recursive = T)
load("./cfrna/results/wgcna_main/cfwgcna.RData")
load("./scrna/data/input/sc_data.RData")
for (i in unique(metacell_seurat$cell_type)) {
#get cell specific data
metacell_subset<-subset(x = metacell_seurat, subset = cell_type %in% i)
#keep 5000 most variable genes
metacell_subset <- FindVariableFeatures(metacell_subset, nfeatures = 5000)
genes.keep <- VariableFeatures(metacell_subset)
metacell_input <- t(as.data.frame(GetAssayData(metacell_subset, assay='RNA', slot='data')[genes.keep,]))
setLabels <- c("cfRNA", "scRNA") 
multiExpr <- list(cfRNA = list(data = cfrna_wgcna_input), scRNA = list(data = metacell_input))
multiColor <- list(cfRNA = cfmergedcolors)
mp_cf_sc <- modulePreservation(multiExpr, multiColor, referenceNetworks = 1, nPermutations = 100, randomSeed = 1, quickCor = 0, maxModuleSize=5000, maxGoldModuleSize=5000, parallelCalculation=T)
ref <- 1 
test <- 2 
statsZ_all <- cbind(mp_cf_sc$quality$Z[[ref]][[test]][, -1], mp_cf_sc$preservation$Z[[ref]][[test]][, -1])
save(statsZ_all, file=paste0("./scrna/results/cf_mp/", i, "_", "statsZ.RData"))
}
```

```
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
```

# Module presevation calculation (exoRNA modules => cell specific scRNA data)
## warning - resource intensive
***


```r
enableWGCNAThreads(nThreads=snakemake@threads) 
```

```
## Allowing parallel execution with up to 15 working processes.
```

```r
dir.create("./scrna/results/exo_mp/", recursive = T)
load("./exorna/results/wgcna_main/exowgcna.RData")
#load("./scrna/data/input/sc_plots.RData")
for (i in unique(metacell_seurat$cell_type)) {
#get cell specific data
metacell_subset<-subset(x = metacell_seurat, subset = cell_type %in% i)
#keep 5000 most variable genes
metacell_subset <- FindVariableFeatures(metacell_subset, nfeatures = 5000)
genes.keep <- VariableFeatures(metacell_subset)
metacell_input <- t(as.data.frame(GetAssayData(metacell_subset, assay='RNA', slot='data')[genes.keep,]))
setLabels <- c("cfRNA", "scRNA") 
multiExpr <- list(exoRNA = list(data = exorna_wgcna_input), scRNA = list(data = metacell_input))
multiColor <- list(exoRNA = exomergedcolors)
mp_exo_sc <- modulePreservation(multiExpr, multiColor, referenceNetworks = 1, nPermutations = 100, randomSeed = 1, quickCor = 0, maxModuleSize=5500, maxGoldModuleSize=5500, parallelCalculation=T)
ref <- 1 
test <- 2 
statsZ_all <- cbind(mp_exo_sc$quality$Z[[ref]][[test]][, -1], mp_exo_sc$preservation$Z[[ref]][[test]][, -1])
save(statsZ_all, file=paste0("./scrna/results/exo_mp/", i, "_", "statsZ.RData"))
}
```

```
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
##   ..checking data for excessive amounts of missing data..
##   ..calculating observed preservation values
##   ..calculating permutation Z scores
##  ..Working with set 1 as reference set
```


```r
sessionInfo()
```

```
## R version 4.1.2 (2021-11-01)
## Platform: x86_64-conda-linux-gnu (64-bit)
## Running under: Debian GNU/Linux 9 (stretch)
## 
## Matrix products: default
## BLAS/LAPACK: /home/si48met/miniconda3/envs/cfwgcna_hcc/lib/libopenblasp-r0.3.21.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] stringr_1.4.1         tidyr_1.2.1           ggrepel_0.9.1        
##  [4] tibble_3.1.8          RColorBrewer_1.1-3    ggplot2_3.3.6        
##  [7] cowplot_1.1.1         WGCNA_1.71            fastcluster_1.2.3    
## [10] dynamicTreeCut_1.63-1 dplyr_1.0.10          forcats_0.5.2        
## [13] sp_1.5-0              SeuratObject_4.1.1    Seurat_4.1.0         
## 
## loaded via a namespace (and not attached):
##   [1] backports_1.4.1        Hmisc_4.7-1            plyr_1.8.7            
##   [4] igraph_1.3.0           lazyeval_0.2.2         splines_4.1.2         
##   [7] listenv_0.8.0          scattermore_0.8        GenomeInfoDb_1.30.0   
##  [10] digest_0.6.29          foreach_1.5.2          htmltools_0.5.3       
##  [13] GO.db_3.14.0           fansi_1.0.3            checkmate_2.1.0       
##  [16] magrittr_2.0.3         memoise_2.0.1          tensor_1.5            
##  [19] cluster_2.1.3          doParallel_1.0.17      ROCR_1.0-11           
##  [22] Biostrings_2.62.0      globals_0.16.1         matrixStats_0.62.0    
##  [25] spatstat.sparse_2.1-1  jpeg_0.1-9             colorspace_2.0-3      
##  [28] blob_1.2.3             xfun_0.32              RCurl_1.98-1.8        
##  [31] crayon_1.5.1           jsonlite_1.8.0         progressr_0.11.0      
##  [34] spatstat.data_2.2-0    impute_1.68.0          survival_3.4-0        
##  [37] zoo_1.8-10             iterators_1.0.14       glue_1.6.2            
##  [40] polyclip_1.10-0        gtable_0.3.1           zlibbioc_1.40.0       
##  [43] XVector_0.34.0         leiden_0.4.2           future.apply_1.9.1    
##  [46] BiocGenerics_0.40.0    abind_1.4-5            scales_1.2.1          
##  [49] DBI_1.1.3              spatstat.random_2.2-0  miniUI_0.1.1.1        
##  [52] Rcpp_1.0.9             htmlTable_2.4.1        viridisLite_0.4.1     
##  [55] xtable_1.8-4           reticulate_1.26        spatstat.core_2.4-4   
##  [58] foreign_0.8-82         bit_4.0.4              preprocessCore_1.56.0 
##  [61] Formula_1.2-4          stats4_4.1.2           htmlwidgets_1.5.4     
##  [64] httr_1.4.4             ellipsis_0.3.2         ica_1.0-3             
##  [67] pkgconfig_2.0.3        nnet_7.3-17            sass_0.4.2            
##  [70] uwot_0.1.11            deldir_1.0-6           utf8_1.2.2            
##  [73] tidyselect_1.1.2       rlang_1.0.5            reshape2_1.4.4        
##  [76] later_1.2.0            AnnotationDbi_1.56.1   munsell_0.5.0         
##  [79] tools_4.1.2            cachem_1.0.6           cli_3.4.0             
##  [82] generics_0.1.3         RSQLite_2.2.8          ggridges_0.5.3        
##  [85] evaluate_0.16          fastmap_1.1.0          yaml_2.3.5            
##  [88] goftest_1.2-3          knitr_1.40             bit64_4.0.5           
##  [91] fitdistrplus_1.1-8     purrr_0.3.4            RANN_2.6.1            
##  [94] KEGGREST_1.34.0        pbapply_1.5-0          future_1.28.0         
##  [97] nlme_3.1-159           mime_0.12              rstudioapi_0.14       
## [100] compiler_4.1.2         plotly_4.10.0          png_0.1-7             
## [103] spatstat.utils_2.3-1   bslib_0.4.0            stringi_1.7.6         
## [106] rgeos_0.5-9            lattice_0.20-45        Matrix_1.4-1          
## [109] vctrs_0.4.1            pillar_1.8.1           lifecycle_1.0.1       
## [112] spatstat.geom_2.4-0    lmtest_0.9-40          jquerylib_0.1.4       
## [115] RcppAnnoy_0.0.19       bitops_1.0-7           data.table_1.14.2     
## [118] irlba_2.3.5            httpuv_1.6.6           patchwork_1.1.2       
## [121] latticeExtra_0.6-30    R6_2.5.1               promises_1.2.0.1      
## [124] KernSmooth_2.23-20     gridExtra_2.3          IRanges_2.28.0        
## [127] parallelly_1.32.1      codetools_0.2-18       MASS_7.3-58.1         
## [130] assertthat_0.2.1       withr_2.5.0            sctransform_0.3.4     
## [133] GenomeInfoDbData_1.2.7 S4Vectors_0.32.3       mgcv_1.8-40           
## [136] parallel_4.1.2         grid_4.1.2             rpart_4.1.16          
## [139] rmarkdown_2.16         Rtsne_0.16             base64enc_0.1-3       
## [142] Biobase_2.54.0         shiny_1.7.2            interp_1.1-3
```
