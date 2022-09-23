---
title: "WGCNA of cfRNA samples & generation of Suppl. File 1,3"
subtitle: "Data of Zhu et al. https://doi.org/10.7150%2Fthno.48206"
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
library(igraph)
library(RColorBrewer)
library(tibble)
library(ggrepel)
library(tidyr)
library(ragg)
library(cowplot)
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(openxlsx)
library(stringr)
```

## Data cleanup 
***


```r
#load gene expression matrix and metadata
load("./cfrna/data/input/cfrna_input.RData")
#check with WGCNA function if there are any bad samples or genes (e.g. too many missing values)
gsg <- goodSamplesGenes(cfrna_wgcna)
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
sampleTree <- hclust(dist(cfrna_wgcna), method = "average")
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
##  0  1  2 
##  6 32 27
```

```r
#choose only non-outlier samples
cfrna_wgcna_input<-cfrna_wgcna[clust %in% c(1,2), ]
#check data dimensions 
dim(cfrna_wgcna_input)
```

```
## [1]   59 9037
```

# WGCNA 
## picking a power value
***


```r
#use multiple threads 
enableWGCNAThreads(nThreads=snakemake@threads) 
```

```
## Allowing parallel execution with up to 15 working processes.
```

```r
#picking soft threshold power
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
#calculate the fit to a scale-free network for different powers 
sft <- pickSoftThreshold(cfrna_wgcna_input, powerVector = powers, networkType = "unsigned")
```

```
##    Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
## 1      1   0.6790  1.4100          0.947  3040.0  3040.000   4690
## 2      2   0.0199 -0.0647          0.746  1470.0  1340.000   3140
## 3      3   0.8100 -0.5660          0.910   845.0   663.000   2320
## 4      4   0.9560 -0.7970          0.970   536.0   355.000   1810
## 5      5   0.9750 -0.9160          0.975   363.0   202.000   1470
## 6      6   0.9840 -0.9760          0.981   259.0   120.000   1220
## 7      7   0.9800 -1.0100          0.975   192.0    73.600   1040
## 8      8   0.9790 -1.0300          0.974   146.0    46.400    896
## 9      9   0.9740 -1.0500          0.968   114.0    29.800    784
## 10    10   0.9740 -1.0600          0.967    91.3    19.700    694
## 11    12   0.9730 -1.0600          0.968    61.1     8.910    557
## 12    14   0.9640 -1.0700          0.958    43.1     4.280    459
## 13    16   0.9590 -1.0800          0.956    31.6     2.150    385
## 14    18   0.9490 -1.0900          0.947    23.9     1.120    332
## 15    20   0.9410 -1.1000          0.946    18.6     0.616    289
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
  ggtitle("Scale independence in cfRNA WGCNA")
#plot 2; mean connectivity (the mean of connections of genes in the network) generally should be about or below lower hundreds
conne<-ggplot(sft$fitIndices, aes(x=Power, y=mean.k., label=Power)) +
  geom_point(size=8) +
  scale_y_log10() + 
  xlab("Soft Threshold (power)") +
  ggtitle("Mean Connectivity") + 
  theme_minimal(base_size = 26, base_family = "Arial") +
  theme(plot.title = element_text(color="black", size=28, face="bold", hjust = 0.5)) +
  ggtitle("Mean connectivity in cfRNA WGCNA")
cfrna_param<-ggdraw() + draw_plot(rsq, x=0, y=0, width = .5, height=1) + draw_plot(conne, x=0.5, y=0, width = .5, height=1)
plot(cfrna_param)
```

![](markdown_results/power-1.png)<!-- -->

## main steps of WGCNA
* build adjacency matrix with the power value 
* calculate TOM and transform it to a disTOM
* hierarchical clustering
* merging similar modules
***


```r
#based on our results and WGCNS recommendations we choose a power value of 6
softPower<-6
enableWGCNAThreads(nThreads=snakemake@threads)
```

```
## Allowing parallel execution with up to 15 working processes.
```

```r
#building an adjacency matrix
adjacency <- adjacency(cfrna_wgcna_input, power = softPower, type = "unsigned")
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
##  ..cutHeight not given, setting it to 0.995  ===>  99% of the (truncated) height range in dendro.
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
##     black      blue     brown     green      grey   magenta      pink    purple 
##       261      1640       922       338        45       159       238       158 
##       red turquoise    yellow 
##       316      4620       340
```

```r
#calculate the module eigengenes (1st principal component of modules)
MEList <- moduleEigengenes(cfrna_wgcna_input, colors = dynamicColors)
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

![](markdown_results/modules-1.png)<!-- -->

```r
#call an automatic merging function
merge <- mergeCloseModules(cfrna_wgcna_input, dynamicColors, cutHeight = MEDissThres)
```

```
##  mergeCloseModules: Merging modules whose distance is less than 0.2
##    Calculating new MEs...
```

```r
#merged module colors
cfmergedcolors <- merge$colors
# get the eigengenes of the new merged modules
mergedMEs <- merge$newMEs
#plot the gene dendrogram to compare the old and new modules
plotDendroAndColors(geneTree, cbind(dynamicColors, cfmergedcolors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
```

![](markdown_results/modules-2.png)<!-- -->

```r
n<-0
for (i in unique(cfmergedcolors)){
n<-n+1
if (i=="grey") {
} else{
cfmergedcolors[cfmergedcolors==i]<-paste0("cf-", "", i)
} 
}
#get a list of new modules 
table(cfmergedcolors)
```

```
## cfmergedcolors
##     cf-black      cf-blue     cf-brown   cf-magenta    cf-purple       cf-red 
##          599         1640         1160          159          158          316 
## cf-turquoise    cf-yellow         grey 
##         4620          340           45
```

```r
#plot the new merged modules by size
df<-data.frame(dc=cfmergedcolors) %>% 
  dplyr::count(dc) %>% mutate(col=str_remove(dc, "cf-"))
cfrna_bar<-ggplot(df, aes(x=reorder(dc, n), y=n, fill=dc)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = as.character(df$col)) +
  theme_minimal(base_size = 26, base_family = "Arial") +
  theme(legend.position = 'None', panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), axis.title=element_blank(), plot.title = element_text(color="black", size=28, face="bold", hjust = 0.5)) +
  ggtitle("Number of genes in cfRNA modules") +
  scale_y_continuous(expand = c(0, 0))
plot(cfrna_bar)
```

![](markdown_results/modules-3.png)<!-- -->

```r
#get the top hub gene of each module (most connected gene in a module)
chooseTopHubInEachModule(cfrna_wgcna_input, cfmergedcolors, power=softPower, type="unsigned")
```

```
##     cf-black      cf-blue     cf-brown   cf-magenta    cf-purple       cf-red 
##     "RHBDF2"       "RPL5"      "AGAP3"     "CACNB1"     "STK17B"      "KDM5A" 
## cf-turquoise    cf-yellow 
##       "TLK1"       "TCF4"
```

```r
#filter metadata appropriately (no outliers)
cfrna_meta_filt<-cfrna_metadata[clust %in% c(1,2), ]
#clean the metadata (for downstream correlation to be possible) and choose the columns of interest to correlate with modules 
#gender assignments are transformed from {male/female} => {0/1}
cfrna_meta_filt$gender<-ifelse(cfrna_meta_filt$gender=="male", 0, 1)
#disease state assignments are transformed from {cancer/healthy} => {1/0}
cfrna_meta_filt$disease_state<-ifelse(cfrna_meta_filt$disease_state=="liver cancer patient", 1, 0)
#*choose columns of interest
cfrna_meta_input<- cfrna_meta_filt %>% 
  dplyr::select(age, gender, disease_state)
dir.create("./cfrna/results/")
dir.create("./cfrna/results/wgcna_main/")
save(cfrna_wgcna_input, cfrna_meta_input, cfmergedcolors, cfrna_param, cfrna_bar, file="./cfrna/results/wgcna_main/cfwgcna.RData")
```

# Exporting modules (Cytoscape preparation)
***


```r
dir.create("./cfrna/results/modules/")
options(stringsAsFactors = FALSE)
probes <- colnames(cfrna_wgcna_input)
for (modules in c("cf-blue", "cf-turquoise")) {
  inModule <- is.finite(match(cfmergedcolors, modules))
  modProbes <- probes[inModule]
  nTop <- 30 
  IMConn <- softConnectivity(cfrna_wgcna_input[, modProbes], type = "unsigned", power = 6)
  top <- (rank(-IMConn) <= nTop)
  modTOM <- TOM[inModule, inModule]
  cyt <- exportNetworkToCytoscape(modTOM[top, top],
                                 edgeFile = paste("./cfrna/results/modules/edges_", modules, ".txt", sep=""),
                                 nodeFile = paste("./cfrna/results/modules/nodes_", modules, ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.1,
                                 nodeNames = modProbes[top],
                                 nodeAttr = cfmergedcolors[inModule][top])
}
```

```
##  softConnectivity: FYI: connecitivty of genes with less than 20 valid samples will be returned as NA.
##  ..calculating connectivities.. 
##  softConnectivity: FYI: connecitivty of genes with less than 20 valid samples will be returned as NA.
##  ..calculating connectivities..
```

# Generation of Supplementary File 1
***


```r
#export 50 top hub genes per module based on absolute values of gene connectivity 
dataset_names <- list()
m <- 0
module <- unique(cfmergedcolors)
module <- module[module != "grey"]
for (modules in module) {
  m <- m + 1
  inModule <- is.finite(match(cfmergedcolors, modules))
  modProbes <- probes[inModule]
  IMConn <- softConnectivity(cfrna_wgcna_input[, modProbes], type = "unsigned", power = 6)
  df_temp <- data.frame(genes=modProbes, module = modules, connectivity=IMConn)
  dataset_names[[m]] <- df_temp %>% top_n(., 50, connectivity) %>% arrange(., desc(connectivity))
}
```

```
##  softConnectivity: FYI: connecitivty of genes with less than 20 valid samples will be returned as NA.
##  ..calculating connectivities.. 
##  softConnectivity: FYI: connecitivty of genes with less than 20 valid samples will be returned as NA.
##  ..calculating connectivities.. 
##  softConnectivity: FYI: connecitivty of genes with less than 20 valid samples will be returned as NA.
##  ..calculating connectivities.. 
##  softConnectivity: FYI: connecitivty of genes with less than 20 valid samples will be returned as NA.
##  ..calculating connectivities.. 
##  softConnectivity: FYI: connecitivty of genes with less than 20 valid samples will be returned as NA.
##  ..calculating connectivities.. 
##  softConnectivity: FYI: connecitivty of genes with less than 20 valid samples will be returned as NA.
##  ..calculating connectivities.. 
##  softConnectivity: FYI: connecitivty of genes with less than 20 valid samples will be returned as NA.
##  ..calculating connectivities.. 
##  softConnectivity: FYI: connecitivty of genes with less than 20 valid samples will be returned as NA.
##  ..calculating connectivities..
```

```r
dir.create("./suppl_files/")
write.xlsx(dataset_names, sheetName = module, file = './suppl_files/supplementary_file1.xlsx') 
```

# Generation of Supplementary File 3
***


```r
#GO Biological process enrichment analysis per module
dataset_path <- list()
x <- 0
module_genes<-data.frame(module = cfmergedcolors, gene = colnames(cfrna_wgcna_input)) 
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
write.xlsx(dataset_path, sheetName = module, file = './suppl_files/supplementary_file3.xlsx') 
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
##  [1] stringr_1.4.1         openxlsx_4.2.5        clusterProfiler_4.2.0
##  [4] org.Hs.eg.db_3.14.0   AnnotationDbi_1.56.1  IRanges_2.28.0       
##  [7] S4Vectors_0.32.3      Biobase_2.54.0        BiocGenerics_0.40.0  
## [10] ragg_1.2.2            tidyr_1.2.1           ggrepel_0.9.1        
## [13] tibble_3.1.8          RColorBrewer_1.1-3    igraph_1.3.0         
## [16] dendextend_1.16.0     ggplot2_3.3.6         cowplot_1.1.1        
## [19] WGCNA_1.71            fastcluster_1.2.3     dynamicTreeCut_1.63-1
## [22] dplyr_1.0.10          forcats_0.5.2        
## 
## loaded via a namespace (and not attached):
##   [1] shadowtext_0.1.2       backports_1.4.1        Hmisc_4.7-1           
##   [4] fastmatch_1.1-3        systemfonts_1.0.4      plyr_1.8.7            
##   [7] lazyeval_0.2.2         splines_4.1.2          BiocParallel_1.28.3   
##  [10] GenomeInfoDb_1.30.0    digest_0.6.29          foreach_1.5.2         
##  [13] yulab.utils_0.0.5      htmltools_0.5.3        GOSemSim_2.20.0       
##  [16] viridis_0.6.2          GO.db_3.14.0           fansi_1.0.3           
##  [19] magrittr_2.0.3         checkmate_2.1.0        memoise_2.0.1         
##  [22] cluster_2.1.3          doParallel_1.0.17      Biostrings_2.62.0     
##  [25] graphlayouts_0.8.1     matrixStats_0.62.0     enrichplot_1.14.1     
##  [28] jpeg_0.1-9             colorspace_2.0-3       blob_1.2.3            
##  [31] textshaping_0.3.6      xfun_0.32              crayon_1.5.1          
##  [34] RCurl_1.98-1.8         jsonlite_1.8.0         scatterpie_0.1.8      
##  [37] impute_1.68.0          survival_3.4-0         iterators_1.0.14      
##  [40] ape_5.6-2              glue_1.6.2             polyclip_1.10-0       
##  [43] gtable_0.3.1           zlibbioc_1.40.0        XVector_0.34.0        
##  [46] scales_1.2.1           DOSE_3.20.0            DBI_1.1.3             
##  [49] Rcpp_1.0.9             viridisLite_0.4.1      htmlTable_2.4.1       
##  [52] tidytree_0.4.0         gridGraphics_0.5-1     foreign_0.8-82        
##  [55] bit_4.0.4              preprocessCore_1.56.0  Formula_1.2-4         
##  [58] htmlwidgets_1.5.4      httr_1.4.4             fgsea_1.20.0          
##  [61] ellipsis_0.3.2         pkgconfig_2.0.3        farver_2.1.1          
##  [64] nnet_7.3-17            sass_0.4.2             deldir_1.0-6          
##  [67] utf8_1.2.2             labeling_0.4.2         ggplotify_0.1.0       
##  [70] tidyselect_1.1.2       rlang_1.0.5            reshape2_1.4.4        
##  [73] munsell_0.5.0          tools_4.1.2            cachem_1.0.6          
##  [76] downloader_0.4         cli_3.4.0              generics_0.1.3        
##  [79] RSQLite_2.2.8          evaluate_0.16          fastmap_1.1.0         
##  [82] yaml_2.3.5             ggtree_3.2.0           knitr_1.40            
##  [85] bit64_4.0.5            tidygraph_1.2.2        zip_2.2.1             
##  [88] purrr_0.3.4            KEGGREST_1.34.0        ggraph_2.0.6          
##  [91] nlme_3.1-159           aplot_0.1.7            DO.db_2.9             
##  [94] compiler_4.1.2         rstudioapi_0.14        png_0.1-7             
##  [97] treeio_1.18.0          tweenr_2.0.2           bslib_0.4.0           
## [100] stringi_1.7.6          highr_0.9              lattice_0.20-45       
## [103] Matrix_1.4-1           vctrs_0.4.1            pillar_1.8.1          
## [106] lifecycle_1.0.1        jquerylib_0.1.4        data.table_1.14.2     
## [109] bitops_1.0-7           patchwork_1.1.2        qvalue_2.26.0         
## [112] R6_2.5.1               latticeExtra_0.6-30    gridExtra_2.3         
## [115] codetools_0.2-18       MASS_7.3-58.1          assertthat_0.2.1      
## [118] withr_2.5.0            GenomeInfoDbData_1.2.7 parallel_4.1.2        
## [121] grid_4.1.2             rpart_4.1.16           ggfun_0.0.7           
## [124] rmarkdown_2.16         ggforce_0.3.4          base64enc_0.1-3       
## [127] interp_1.1-3
```
