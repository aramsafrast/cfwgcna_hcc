---
title: "WGCNA of exoRNA samples & generation of Suppl. File 2,4"
subtitle: "Data of Li et al. https://doi.org/10.1093/nar/gkx891"
author: "Aram Safrastyan"
date: "19 September, 2022"
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
library(forcats)
library(dplyr)
library(WGCNA)
library(cowplot)
library(ggplot2)
library(dendextend)
library(RColorBrewer)
library(tibble)
library(ggrepel)
library(tidyr)
library(ragg)
library(cowplot)
library(stringr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(openxlsx)
```

# Data cleanup 
***


```r
#load gene expression matrix and metadata
load("./exorna/data/input/exorna_input.RData")
#check with WGCNA function if there are any bad samples or genes (e.g. too many missing values)
gsg <- goodSamplesGenes(exorna_wgcna)
```

```
##  Flagging genes and samples with too many missing values...
##   ..step 1
```

```r
gsg$allOK #no bad samples of genes detected
```

```
## [1] TRUE
```

```r
#for outlier detection hierarchical clustering is done
sampleTree <- hclust(dist(exorna_wgcna), method = "average")
#plot the results as a dendrogram 
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
```

![](markdown_results/clean-1.png)<!-- -->

```r
#based on dendrogram choose cutting height - the outlier samples to filter out 
clust <- cutreeStatic(sampleTree, cutHeight = 100, minSize = 2)
#examine the clusters - one should be very small and contain only the outlier samples
table(clust)
```

```
## clust
##  1 
## 21
```

```r
#choose only non-outlier samples
exorna_wgcna_input<-exorna_wgcna[clust==1, ]
#check data dimensions 
dim(exorna_wgcna_input)
```

```
## [1]    21 11481
```

# WGCNA 
## picking a power value
***


```r
#use 8 threads 
enableWGCNAThreads(nThreads=snakemake@threads) 
```

```
## Allowing parallel execution with up to 15 working processes.
```

```r
#picking soft threshold power
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
#calculate the fit to a scale-free network for different powers 
sft <- pickSoftThreshold(exorna_wgcna_input, powerVector = powers, networkType = "unsigned")
```

```
##    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
## 1      1   0.0174 -0.374          0.361 2690.00  2.51e+03   4390
## 2      2   0.7290 -1.330          0.745  981.00  8.02e+02   2440
## 3      3   0.9580 -1.480          0.946  451.00  3.08e+02   1610
## 4      4   0.9630 -1.460          0.956  241.00  1.34e+02   1190
## 5      5   0.9480 -1.410          0.936  145.00  6.37e+01    939
## 6      6   0.9430 -1.350          0.929   94.20  3.25e+01    772
## 7      7   0.9240 -1.310          0.902   65.40  1.75e+01    654
## 8      8   0.9220 -1.270          0.901   47.70  9.90e+00    565
## 9      9   0.9190 -1.240          0.900   36.20  5.82e+00    496
## 10    10   0.9200 -1.210          0.906   28.30  3.54e+00    441
## 11    12   0.9300 -1.170          0.925   18.60  1.42e+00    357
## 12    14   0.9250 -1.140          0.933   13.00  6.16e-01    296
## 13    16   0.9220 -1.130          0.937    9.57  2.89e-01    250
## 14    18   0.9320 -1.120          0.954    7.26  1.44e-01    214
## 15    20   0.9290 -1.120          0.959    5.65  7.55e-02    185
```

```r
#plotting the results; plot 1
rsq<-ggplot(sft$fitIndices, aes(x=Power, y=(-sign(slope)*SFT.R.sq), label=Power)) + 
  geom_text(col="blue", size=10, family="Arial", fontface="bold") +
  ylab("Scale Free Topology Model Fit signed Rˆ2") +
  xlab("Soft Threshold (power)") + 
  geom_hline(yintercept = 0.85) + 
  theme_minimal(base_size = 26, base_family = "Arial") +
  theme(plot.title = element_text(color="black", size=28, face="bold", hjust = 0.5)) +
  ggtitle("Scale independence in exoRNA WGCNA")
#plot 2; mean connectivity (the mean of connections of genes in the network) generally should be about or below lower hundreds
conne<-ggplot(sft$fitIndices, aes(x=Power, y=mean.k., label=Power)) +
  geom_point(size=8) +
  scale_y_log10() + 
  xlab("Soft Threshold (power)") +
  ggtitle("Mean Connectivity") + 
  theme_minimal(base_size = 26, base_family = "Arial") +
  theme(plot.title = element_text(color="black", size=28, face="bold", hjust = 0.5)) +
  ggtitle("Mean connectivity in exoRNA WGCNA")
exorna_param<-ggdraw() + draw_plot(rsq, x=0, y=0, width = .5, height=1) + draw_plot(conne, x=0.5, y=0, width = .5, height=1) 
plot(exorna_param)
```

<img src="markdown_results/power-1.png" style="display: block; margin: auto;" />

# main steps of WGCNA
* build adjacency matrix with the power value 
* calculate TOM and transform it to a disTOM
* hierarchical clustering
* merging similar modules
***


```r
#based on our results and WGCNS recommendations we choose a power value of 8
softPower<-8
#building an adjacency matrix
adjacency <- adjacency(exorna_wgcna_input, power = softPower, type = "unsigned")
#turn adjacency into topological overlap matrix (TOM) 
TOM <- TOMsimilarity(adjacency, TOMType = "unsigned")
```

```
## ..connectivity..
## ..matrix multiplication (system BLAS)..
## ..normalization..
## ..done.
```

```r
#transformation into a dissimilarity measure 
dissTOM <- 1-TOM
# Call the hierarchical clustering function and use the dissimilarity measure 
geneTree <- hclust(as.dist(dissTOM), method = "average")
#we pick a minimal threshold of 30 genes per potential module 
minModuleSize <- 30
# Module identification using dynamic tree cut:
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, minClusterSize = minModuleSize, deepSplit=2)
```

```
##  ..cutHeight not given, setting it to 0.996  ===>  99% of the (truncated) height range in dendro.
##  ..done.
```

```r
# Convert numeric labels into colors
dynamicColors <- labels2colors(dynamicMods)
#count the resulting modules
table(dynamicColors)
```

```
## dynamicColors
##      blue     brown     green      grey       red turquoise    yellow 
##       670       411       112      4941        44      5055       248
```

```r
#calculate the module eigengenes (1st principal component of modules)
MEList <- moduleEigengenes(exorna_wgcna_input, colors = dynamicColors)
MEs <- MEList$eigengenes
#calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs)
#cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")
#plot the result
d<-as.dendrogram(METree)
lab<-METree$labels[order.dendrogram(d)]
d<-color_labels(d ,labels=lab, col = sub("ME","",lab) )
par(mar=c(5,4,4,6)) 
plot(d, 
     main = "Clustering of module eigengenes", 
     horiz=T, 
     xlab = "", sub = "",cex=0.5)
#we pick a dissimilarity threshold of 0.2 (=correlation of 0.8) to merge very similar modules
MEDissThres <- 0.2
#plot the threshold line
abline(v=0.2, col=3)
```

<img src="markdown_results/modules-1.png" style="display: block; margin: auto;" />

```r
#call an automatic merging function
merge <- mergeCloseModules(exorna_wgcna_input, dynamicColors, cutHeight = MEDissThres)
```

```
##  mergeCloseModules: Merging modules whose distance is less than 0.2
##    Calculating new MEs...
```

```r
#merged module colors
exomergedcolors <- merge$colors
# get the eigengenes of the new merged modules
mergedMEs <- merge$newMEs
#plot the gene dendrogram to compare the old and new modules
plotDendroAndColors(geneTree, cbind(dynamicColors, exomergedcolors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
```

<img src="markdown_results/modules-2.png" style="display: block; margin: auto;" />

```r
n<-0
for (i in unique(exomergedcolors)){
n<-n+1
if (i=="grey") {
} else{
exomergedcolors[exomergedcolors==i]<-paste0("exo-", "", i)
} 
}
#get a list of new modules 
table(exomergedcolors)
```

```
## exomergedcolors
##      exo-blue     exo-brown     exo-green       exo-red exo-turquoise 
##           670           411           112            44          5055 
##    exo-yellow          grey 
##           248          4941
```

```r
#plot the new merged modules by size
df<-data.frame(dc=exomergedcolors) %>% 
  dplyr::count(dc) %>% mutate(col=str_remove(dc, "exo-"))
exorna_bar<-ggplot(df, aes(x=reorder(dc, n), y=n, fill=dc)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = as.character(df$col)) +
  theme_minimal(base_size = 26, base_family = "Arial") +
  theme(legend.position = 'None', panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), axis.title=element_blank(), plot.title = element_text(color="black", size=28, face="bold", hjust = 0.5)) +
  ggtitle("Number of genes in exoRNA modules") +
  scale_y_continuous(expand = c(0, 0))
plot(exorna_bar)
```

<img src="markdown_results/modules-3.png" style="display: block; margin: auto;" />

```r
#get the top hub gene of each module (most connected gene in a module)
chooseTopHubInEachModule(exorna_wgcna_input, exomergedcolors, power=softPower, type="unsigned")
```

```
##      exo-blue     exo-brown     exo-green       exo-red exo-turquoise 
##      "DNAH14"       "RPL34"      "DCAF12"      "BNIP3L"        "MTPN" 
##    exo-yellow 
##       "AKAP7"
```

```r
#filter metadata appropriately (no outliers)
exorna_meta_input<-exorna_metadata[clust==1, ]
dir.create("./exorna/results/")
dir.create("./exorna/results/wgcna_main/")
save(exorna_wgcna_input, exorna_meta_input, exomergedcolors, exorna_param, exorna_bar, file="./exorna/results/wgcna_main/exowgcna.RData")
```

# Generation of Supplementary File 2
***


```r
#export 50 top hub genes per module based on absolute values of module membership
dataset_names <- list()
m <- 0
probes <- colnames(exorna_wgcna_input)
module <- unique(exomergedcolors)
module <- module[module != "grey"]
for (modules in module) {
  m <- m + 1
  inModule <- is.finite(match(exomergedcolors, modules))
  modProbes <- probes[inModule]
  IMConn <- softConnectivity(exorna_wgcna_input[, modProbes], type = "unsigned", power = 8)
  df_temp <- data.frame(genes=modProbes, module = modules, connectivity=IMConn)
  dataset_names[[m]] <- df_temp %>% top_n(., 50, connectivity) %>% arrange(., desc(connectivity))
}
```

```
##  softConnectivity: FYI: connecitivty of genes with less than 7 valid samples will be returned as NA.
##  ..calculating connectivities.. 
##  softConnectivity: FYI: connecitivty of genes with less than 7 valid samples will be returned as NA.
##  ..calculating connectivities.. 
##  softConnectivity: FYI: connecitivty of genes with less than 7 valid samples will be returned as NA.
##  ..calculating connectivities.. 
##  softConnectivity: FYI: connecitivty of genes with less than 7 valid samples will be returned as NA.
##  ..calculating connectivities.. 
##  softConnectivity: FYI: connecitivty of genes with less than 7 valid samples will be returned as NA.
##  ..calculating connectivities.. 
##  softConnectivity: FYI: connecitivty of genes with less than 7 valid samples will be returned as NA.
##  ..calculating connectivities..
```

```r
dir.create("./suppl_files/")
write.xlsx(dataset_names, sheetName = module, file = './suppl_files/supplementary_file2.xlsx') 
```

# Generation of Supplementary File 4
***


```r
#GO Biological process enrichment analysis per module
dataset_path <- list()
x <- 0
module_genes<-data.frame(module = exomergedcolors, gene = colnames(exorna_wgcna_input)) 
module_genes_entrez<-AnnotationDbi::select(org.Hs.eg.db,
                              keys = module_genes$gene,
                              keytype = "SYMBOL",
                              columns = c("SYMBOL","ENTREZID")) %>% dplyr::left_join(module_genes,by=c("SYMBOL"="gene"))
for (modules in module) {
  x <- x + 1
  mod_genes <- module_genes_entrez %>% dplyr::filter(module==modules) %>% na.omit() %>%
    dplyr::select(ENTREZID) %>% unique()
  go <- enrichGO(gene = mod_genes$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                readable = TRUE)
  go_df<-fortify(go, showCategory = 5) %>% 
    dplyr::select(., c(ID, Description, p.adjust, geneID, Count))
  rownames(go_df) <- NULL
  dataset_path[[x]] <- go_df
}
#export the top 10 results
write.xlsx(dataset_path, sheetName = module, file = './suppl_files/supplementary_file4.xlsx') 
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
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] openxlsx_4.2.5        clusterProfiler_4.2.0 org.Hs.eg.db_3.14.0  
##  [4] AnnotationDbi_1.56.1  IRanges_2.28.0        S4Vectors_0.32.3     
##  [7] Biobase_2.54.0        BiocGenerics_0.40.0   stringr_1.4.1        
## [10] ragg_1.2.2            tidyr_1.2.1           ggrepel_0.9.1        
## [13] tibble_3.1.8          RColorBrewer_1.1-3    dendextend_1.16.0    
## [16] ggplot2_3.3.6         cowplot_1.1.1         WGCNA_1.71           
## [19] fastcluster_1.2.3     dynamicTreeCut_1.63-1 dplyr_1.0.10         
## [22] forcats_0.5.2        
## 
## loaded via a namespace (and not attached):
##   [1] shadowtext_0.1.2       backports_1.4.1        Hmisc_4.7-1           
##   [4] fastmatch_1.1-3        systemfonts_1.0.4      plyr_1.8.7            
##   [7] igraph_1.3.0           lazyeval_0.2.2         splines_4.1.2         
##  [10] BiocParallel_1.28.3    GenomeInfoDb_1.30.0    digest_0.6.29         
##  [13] foreach_1.5.2          yulab.utils_0.0.5      htmltools_0.5.3       
##  [16] GOSemSim_2.20.0        viridis_0.6.2          GO.db_3.14.0          
##  [19] fansi_1.0.3            magrittr_2.0.3         checkmate_2.1.0       
##  [22] memoise_2.0.1          cluster_2.1.3          doParallel_1.0.17     
##  [25] Biostrings_2.62.0      graphlayouts_0.8.1     matrixStats_0.62.0    
##  [28] enrichplot_1.14.1      jpeg_0.1-9             colorspace_2.0-3      
##  [31] blob_1.2.3             textshaping_0.3.6      xfun_0.32             
##  [34] crayon_1.5.1           RCurl_1.98-1.8         jsonlite_1.8.0        
##  [37] scatterpie_0.1.8       impute_1.68.0          survival_3.4-0        
##  [40] iterators_1.0.14       ape_5.6-2              glue_1.6.2            
##  [43] polyclip_1.10-0        gtable_0.3.1           zlibbioc_1.40.0       
##  [46] XVector_0.34.0         scales_1.2.1           DOSE_3.20.0           
##  [49] DBI_1.1.3              Rcpp_1.0.9             viridisLite_0.4.1     
##  [52] htmlTable_2.4.1        tidytree_0.4.0         gridGraphics_0.5-1    
##  [55] foreign_0.8-82         bit_4.0.4              preprocessCore_1.56.0 
##  [58] Formula_1.2-4          htmlwidgets_1.5.4      httr_1.4.4            
##  [61] fgsea_1.20.0           ellipsis_0.3.2         pkgconfig_2.0.3       
##  [64] farver_2.1.1           nnet_7.3-17            sass_0.4.2            
##  [67] deldir_1.0-6           utf8_1.2.2             labeling_0.4.2        
##  [70] ggplotify_0.1.0        tidyselect_1.1.2       rlang_1.0.5           
##  [73] reshape2_1.4.4         munsell_0.5.0          tools_4.1.2           
##  [76] cachem_1.0.6           downloader_0.4         cli_3.4.0             
##  [79] generics_0.1.3         RSQLite_2.2.8          evaluate_0.16         
##  [82] fastmap_1.1.0          yaml_2.3.5             ggtree_3.2.0          
##  [85] knitr_1.40             bit64_4.0.5            tidygraph_1.2.2       
##  [88] zip_2.2.1              purrr_0.3.4            KEGGREST_1.34.0       
##  [91] ggraph_2.0.6           nlme_3.1-159           aplot_0.1.7           
##  [94] DO.db_2.9              compiler_4.1.2         rstudioapi_0.14       
##  [97] png_0.1-7              treeio_1.18.0          tweenr_2.0.2          
## [100] bslib_0.4.0            stringi_1.7.6          highr_0.9             
## [103] lattice_0.20-45        Matrix_1.4-1           vctrs_0.4.1           
## [106] pillar_1.8.1           lifecycle_1.0.1        jquerylib_0.1.4       
## [109] data.table_1.14.2      bitops_1.0-7           patchwork_1.1.2       
## [112] qvalue_2.26.0          R6_2.5.1               latticeExtra_0.6-30   
## [115] gridExtra_2.3          codetools_0.2-18       MASS_7.3-58.1         
## [118] assertthat_0.2.1       withr_2.5.0            GenomeInfoDbData_1.2.7
## [121] parallel_4.1.2         grid_4.1.2             rpart_4.1.16          
## [124] ggfun_0.0.7            rmarkdown_2.16         ggforce_0.3.4         
## [127] base64enc_0.1-3        interp_1.1-3
```
