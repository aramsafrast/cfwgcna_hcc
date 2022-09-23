---
title: "Preprocessing of single-cell RNA samples from livers of healthy donors"
subtitle: "Data of MacParland et al. https://doi.org/10.1038/s41467-018-06318-7"
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
#library(scWGCNA)
library(cowplot)
library(ggplot2)
library(dplyr)
library(Seurat)
library(DT)
library(tibble)
library(tidyr)
library(BiocParallel)
```

#  Download the data (GSE115469)
***


```bash
#download the scRNA count matrix
mkdir -p ./scrna/data/raw/
wget -q -O - https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115469/suppl/GSE115469_Data.csv.gz | gunzip -c > ./scrna/data/raw/scrna_countm.csv 
#download the metadata
wget -q -O - https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115469/suppl/GSE115469_CellClusterType.txt.gz | gunzip -c > ./scrna/data/raw/scrna_metadata_init.txt 
```

# Data cleanup; normalization
***


```r
#set seed for reproducible UMAP calculations
set.seed(1)
#load data into Seurat
ref <- read.csv("./scrna/data/raw/scrna_countm.csv")
ref.txt <- read.delim("./scrna/data/raw/scrna_metadata_init.txt")
#merge cell subtypes
ref.txt$CellType<-ifelse(grepl("Hepatocyte",ref.txt$CellType),"Hepatocyte",ref.txt$CellType)
ref.txt$CellType<-ifelse(grepl("T_Cells",ref.txt$CellType),"T_Cells",ref.txt$CellType)
ref.txt$CellType<-ifelse(grepl("LSECs",ref.txt$CellType),"LSECs",ref.txt$CellType)
ref.txt$CellType<-ifelse(grepl("Macrophage",ref.txt$CellType),"Macrophage",ref.txt$CellType)
ref<-ref %>% column_to_rownames(., var = "X")
sobject <- CreateSeuratObject(ref)
sobject$cell_type<-ref.txt$CellType
#normalize and visualize the data
sobject <- NormalizeData(sobject)
sobject <- FindVariableFeatures(sobject)
sobject <- ScaleData(sobject)
sobject <- RunPCA(sobject, features = VariableFeatures(object = sobject))
sobject <- FindNeighbors(sobject)
sobject <- FindClusters(sobject)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 8444
## Number of edges: 294037
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8950
## Number of communities: 17
## Elapsed time: 1 seconds
```

```r
sobject <- RunUMAP(sobject,  features = VariableFeatures(object = sobject))
sobject@meta.data$cell_type <-recode(sobject@meta.data$cell_type, Hepatocyte="Hepatocytes", T_Cells = "T cells", Macrophage="Macrophages", 'NK-like_Cells'="NK-like cells", Erythroid_Cells="Erythroid cells", Mature_B_Cells="Mature B cells", Plasma_Cells="Plasma cells", Portal_endothelial_Cells="Portal endothelial cells", Hepatic_Stellate_Cells="Hepatic stellate cells")
options(ggrepel.max.overlaps = Inf)
raw_sc_plot<-DimPlot(sobject, reduction = "umap", label = TRUE, repel = TRUE, label.size = 10, group.by = "cell_type", combine = TRUE) +  
  labs(title=paste0("UMAP plot of healthy liver single-cell dataset before metacell transformation", "\n", "(MacParland et al.)"))  +
  theme_classic(base_size = 24) +
  theme(plot.title = element_text(color="black", size=26, face="bold", hjust = 0.5), legend.text=element_text(size=20)) + 
  guides(fill="none", colour = guide_legend(override.aes = list(size=10))) + 
  NoLegend()
plot(raw_sc_plot)
```

![](markdown_results/unnamed-chunk-3-1.png)<!-- -->

# Construct metacells
***


```r
sobject@meta.data %>% dplyr::count(cell_type)
```

```
##                   cell_type    n
## 1            Cholangiocytes  119
## 2           Erythroid cells   93
## 3    Hepatic stellate cells   37
## 4               Hepatocytes 3501
## 5                     LSECs  633
## 6               Macrophages 1192
## 7            Mature B cells  129
## 8             NK-like cells  488
## 9              Plasma cells  511
## 10 Portal endothelial cells  211
## 11                  T cells 1530
```

```r
sobject$metacell_group <- as.character(sobject$cell_type)
#use predefined genes for cell cycle scoring
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
sobject <- CellCycleScoring(sobject, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
```

```
## Warning: The following features are not present in the object: PIMREG, JPT1, not
## searching for symbol synonyms
```

```r
sobject@meta.data %>% dplyr::count(cell_type, Phase)
```

```
##                   cell_type Phase    n
## 1            Cholangiocytes    G1   71
## 2            Cholangiocytes   G2M   22
## 3            Cholangiocytes     S   26
## 4           Erythroid cells    G1   33
## 5           Erythroid cells   G2M   48
## 6           Erythroid cells     S   12
## 7    Hepatic stellate cells    G1   25
## 8    Hepatic stellate cells   G2M    6
## 9    Hepatic stellate cells     S    6
## 10              Hepatocytes    G1 2077
## 11              Hepatocytes   G2M  483
## 12              Hepatocytes     S  941
## 13                    LSECs    G1  410
## 14                    LSECs   G2M  100
## 15                    LSECs     S  123
## 16              Macrophages    G1  639
## 17              Macrophages   G2M  308
## 18              Macrophages     S  245
## 19           Mature B cells    G1   45
## 20           Mature B cells   G2M   38
## 21           Mature B cells     S   46
## 22            NK-like cells    G1  205
## 23            NK-like cells   G2M  109
## 24            NK-like cells     S  174
## 25             Plasma cells    G1   97
## 26             Plasma cells   G2M  252
## 27             Plasma cells     S  162
## 28 Portal endothelial cells    G1  148
## 29 Portal endothelial cells   G2M   29
## 30 Portal endothelial cells     S   34
## 31                  T cells    G1  557
## 32                  T cells   G2M  489
## 33                  T cells     S  484
```

```r
DimPlot(sobject, reduction = "umap", label = TRUE, repel = TRUE, group.by = c("cell_type", "Phase"))
```

![](markdown_results/unnamed-chunk-4-1.png)<!-- -->

```r
#remove hepatic stellate cells due to low number
sobject<-subset(x = sobject, subset = cell_type!="Hepatic stellate cells")
#define the function from the package "scWGCNA"/"hdWGCNA" 
construct_metacells <- function (seurat_obj, name = "agg", k = 50, reduction = "umap", 
    assay = "RNA", slot = "data") 
{
    reduced_coordinates <- as.data.frame(seurat_obj@reductions[[reduction]]@cell.embeddings)
    nn_map <- FNN::knn.index(reduced_coordinates, k = (k - 1))
    row.names(nn_map) <- row.names(reduced_coordinates)
    nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
    good_choices <- seq_len(nrow(nn_map))
    choice <- sample(seq_len(length(good_choices)), size = 1, 
        replace = FALSE)
    chosen <- good_choices[choice]
    good_choices <- good_choices[good_choices != good_choices[choice]]
    it <- 0
    k2 <- k * 2
    get_shared <- function(other, this_choice) {
        k2 - length(union(cell_sample[other, ], this_choice))
    }
    while (length(good_choices) > 0 & it < 5000) {
        it <- it + 1
        choice <- sample(seq_len(length(good_choices)), size = 1, 
            replace = FALSE)
        new_chosen <- c(chosen, good_choices[choice])
        good_choices <- good_choices[good_choices != good_choices[choice]]
        cell_sample <- nn_map[new_chosen, ]
        others <- seq_len(nrow(cell_sample) - 1)
        this_choice <- cell_sample[nrow(cell_sample), ]
        shared <- sapply(others, get_shared, this_choice = this_choice)
        if (max(shared) < 0.9 * k) {
            chosen <- new_chosen
        }
    }
    cell_sample <- nn_map[chosen, ]
    combs <- combn(nrow(cell_sample), 2)
    shared <- apply(combs, 2, function(x) {
        k2 - length(unique(as.vector(cell_sample[x, ])))
    })
    message(paste0("Overlap QC metrics:\nCells per bin: ", k, 
        "\nMaximum shared cells bin-bin: ", max(shared), "\nMean shared cells bin-bin: ", 
        mean(shared), "\nMedian shared cells bin-bin: ", median(shared)))
    if (mean(shared)/k > 0.1) 
        warning("On average, more than 10% of cells are shared between paired bins.")
    exprs_old <- GetAssayData(seurat_obj, assay = assay, slot = slot)
    mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(exprs_old)) %in% 
        cell_sample[x, , drop = FALSE])
    mask <- Matrix::Matrix(mask)
    new_exprs <- (exprs_old %*% mask)/k
    colnames(new_exprs) <- paste0(name, "_", 1:ncol(new_exprs))
    rownames(cell_sample) <- paste0(name, "_", 1:ncol(new_exprs))
    colnames(cell_sample) <- paste0("knn_", 1:ncol(cell_sample))
    seurat_aggr <- CreateSeuratObject(counts = new_exprs)
    seurat_aggr
}
#construct metacells
seurat_list <- list()
for(group in unique(sobject$cell_type)){
  print(group)
  cur_seurat <- subset(sobject, cell_type == group)
  #cur_seurat <- cur_seurat[genes.keep,]
  k<-ifelse(ncol(cur_seurat@assays$RNA) <300, 8, 20)
  cur_metacell_seurat <- construct_metacells(
    cur_seurat, name=group,
    k=k, reduction='umap',
    assay='RNA', slot='data'
  )
  cur_metacell_seurat$cell_type <- as.character(unique(cur_seurat$cell_type))
  seurat_list[[group]] <- cur_metacell_seurat
}
```

```
## [1] "LSECs"
```

```
## Overlap QC metrics:
## Cells per bin: 20
## Maximum shared cells bin-bin: 17
## Mean shared cells bin-bin: 0.660030237924723
## Median shared cells bin-bin: 0
```

```
## [1] "Cholangiocytes"
```

```
## Overlap QC metrics:
## Cells per bin: 8
## Maximum shared cells bin-bin: 7
## Mean shared cells bin-bin: 0.544356435643564
## Median shared cells bin-bin: 0
```

```
## [1] "Macrophages"
```

```
## Overlap QC metrics:
## Cells per bin: 20
## Maximum shared cells bin-bin: 17
## Mean shared cells bin-bin: 0.335514192422087
## Median shared cells bin-bin: 0
```

```
## [1] "T cells"
```

```
## Overlap QC metrics:
## Cells per bin: 20
## Maximum shared cells bin-bin: 17
## Mean shared cells bin-bin: 0.260206804154405
## Median shared cells bin-bin: 0
```

```
## [1] "NK-like cells"
```

```
## Overlap QC metrics:
## Cells per bin: 20
## Maximum shared cells bin-bin: 17
## Mean shared cells bin-bin: 0.840881621045148
## Median shared cells bin-bin: 0
```

```
## [1] "Hepatocytes"
```

```
## Overlap QC metrics:
## Cells per bin: 20
## Maximum shared cells bin-bin: 17
## Mean shared cells bin-bin: 0.111099710904164
## Median shared cells bin-bin: 0
```

```
## [1] "Portal endothelial cells"
```

```
## Overlap QC metrics:
## Cells per bin: 8
## Maximum shared cells bin-bin: 7
## Mean shared cells bin-bin: 0.298165428807211
## Median shared cells bin-bin: 0
```

```
## [1] "Mature B cells"
```

```
## Overlap QC metrics:
## Cells per bin: 8
## Maximum shared cells bin-bin: 7
## Mean shared cells bin-bin: 0.488461538461538
## Median shared cells bin-bin: 0
```

```
## [1] "Plasma cells"
```

```
## Overlap QC metrics:
## Cells per bin: 20
## Maximum shared cells bin-bin: 17
## Mean shared cells bin-bin: 0.795235913879982
## Median shared cells bin-bin: 0
```

```
## [1] "Erythroid cells"
```

```
## Overlap QC metrics:
## Cells per bin: 8
## Maximum shared cells bin-bin: 7
## Mean shared cells bin-bin: 0.669215291750503
## Median shared cells bin-bin: 0
```

```r
# merge all of the metacells objects
metacell_seurat <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)])
#size of metacell seurat file
dim(metacell_seurat)
```

```
## [1] 20007  5266
```

```r
#normalize and visualize the new data
metacell_seurat <- NormalizeData(metacell_seurat)
all.genes <- rownames(metacell_seurat)
metacell_seurat <- ScaleData(metacell_seurat, features = all.genes)
```

```
## Centering and scaling data matrix
```

```r
metacell_seurat <- FindVariableFeatures(metacell_seurat)
metacell_seurat <- RunPCA(metacell_seurat, features = VariableFeatures(object = metacell_seurat))
```

```
## PC_ 1 
## Positive:  ARHGDIB, GMFG, FXYD5, SRGN, CYBA, HLA-E, GSTP1, CLIC1, HLA-A, ARPC1B 
## 	   TMSB10, SLC25A6, IQGAP1, SH3BGRL3, CD74, CLEC2B, TMSB4X, VIM, RPS4Y1, EMP3 
## 	   CARD16, PLAC8, GYPC, LAPTM5, JUNB, GIMAP4, HNRNPA1, S100A6, LDHB, COTL1 
## Negative:  SLC22A1, BAAT, CPB2, IGFBP1, AADAC, CYP8B1, APOA5, LEAP2, CYP4A11, HSD11B1 
## 	   PAH, SERPINA4, GJB1, HPN, PROC, UGT2B4, KHK, ITIH1, CYP2D6, CPS1 
## 	   ACSM2B, MYO1B, ITIH2, BHMT, MASP2, C5, AFM, PLG, TAT-AS1, RDH16 
## PC_ 2 
## Positive:  CD48, HCST, RAC2, ACAP1, TRAF3IP3, LSP1, CORO1A, CXCR4, TBC1D10C, PTPRC 
## 	   CD37, CYTIP, IL2RG, LCP1, IKZF1, RUNX3, SASH3, ALOX5AP, ITGAL, DUSP2 
## 	   ITGB2, RP11-347P5.1, BIN2, CD53, ARHGAP30, IKZF3, STK17B, PTPN7, SLC38A1, TAGAP 
## Negative:  MEIS2, MYCT1, EGFL7, NR2F2, RUNX1T1, NPR1, NOSTRIN, TINAGL1, PTPRB, HSPG2 
## 	   MMRN2, COX7A1, TM4SF1, SOX18, WBP5, TSPAN7, EMCN, CALCRL, PCAT19, TIMP3 
## 	   SRPX, PTRF, CRIP2, IL33, ESAM, RAMP2, HYAL2, GJA4, TEK, PDE2A 
## PC_ 3 
## Positive:  MAFB, MPEG1, TGFBI, CD163, C5AR1, CD68, CYBB, CLEC7A, CSF1R, VSIG4 
## 	   LILRB2, IGSF6, LRRC25, CPVL, TLR2, HCK, SPI1, MS4A7, HMOX1, ADAP2 
## 	   CFD, DMXL2, KCNE3, LINC01272, FPR1, LINC01503, SLC7A7, PTAFR, FCGR1A, NCF2 
## Negative:  IFITM1, ITM2A, ETS1, PRKCH, SKAP1, GZMM, CHST12, TRBC2, PTPN4, FKBP11 
## 	   CTSW, LCK, MATK, PLEKHF1, PRF1, PYHIN1, ZAP70, LBH, CD69, CD2 
## 	   CD247, IKZF3, SH2D1A, TRBC1, LINC00861, GZMA, IL32, SAMD3, LAT, CD7 
## PC_ 4 
## Positive:  CLDN10, KRT19, CLDN4, ANXA13, EPCAM, RAB25, TACSTD2, PLPP2, CHST4, EHF 
## 	   FXYD2, MMP7, KRT7, SLC5A1, USH1C, SMIM22, SLC6A19, SOX9, PRR15L, BICC1 
## 	   B3GNT3, SEZ6L2, CA9, FAM150B, AP1M2, KCNJ15, FAM3B, ERICH5, CFTR, SFRP5 
## Negative:  SOCS2, A2M, TIMP3, A1BG, GPX3, NID1, HP, MT2A, IGFBP2, C9 
## 	   FSTL1, HPX, PLXND1, CD4, CYP2E1, TMEM204, SHANK3, SH3BP5, AHSG, GIMAP4 
## 	   GIMAP7, TSPAN7, LIMS2, SCARF1, NNMT, FAM167B, FHL1, ITPRIP, PAM, ETS1 
## PC_ 5 
## Positive:  MYBL2, CPNE5, POU2AF1, TNFRSF13B, RRM2, MZB1, TNFRSF17, DERL3, JCHAIN, PNOC 
## 	   IGLL5, FCRL5, SHCBP1, KIAA0125, RP11-1070N10.3, UBE2T, ABCB9, SLC38A5, IGLV6-57, MKI67 
## 	   CDKN3, IGLC3, TOP2A, HMBS, ZWINT, CD79A, CENPW, RP11-16E12.2, IGHA2, IGKV1-12 
## Negative:  GIMAP4, GIMAP7, GIMAP1, ANXA1, PLEKHF1, SLFN5, SKAP1, MYO1F, PTGER2, IFITM1 
## 	   PTGDR, MATK, DENND2D, FYB, GZMM, KLRD1, SH2D1A, CD160, APBA2, IFNG 
## 	   SPON2, DOK2, CTSW, CD7, TRBC1, PRKCH, FCRL6, ZAP70, CCL3L3, TRGC1
```

```r
metacell_seurat <- FindNeighbors(metacell_seurat)
```

```
## Computing nearest neighbor graph
```

```
## Computing SNN
```

```r
metacell_seurat <- FindClusters(metacell_seurat)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 5266
## Number of edges: 137921
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9334
## Number of communities: 28
## Elapsed time: 0 seconds
```

```r
metacell_seurat <- RunUMAP(metacell_seurat,  features = VariableFeatures(object = metacell_seurat))
```

```
## 01:43:06 UMAP embedding parameters a = 0.9922 b = 1.112
```

```
## 01:43:07 Read 5266 rows and found 2000 numeric columns
```

```
## 01:43:07 Using Annoy for neighbor search, n_neighbors = 30
```

```
## 01:43:07 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 01:43:11 Writing NN index file to temp file /tmp/RtmpRUEPKg/fileb2f4594eab67
## 01:43:11 Searching Annoy index using 1 thread, search_k = 3000
## 01:43:57 Annoy recall = 100%
## 01:43:58 Commencing smooth kNN distance calibration using 1 thread
## 01:43:59 Found 2 connected components, falling back to 'spca' initialization with init_sdev = 1
## 01:43:59 Initializing from PCA
## 01:43:59 Using 'irlba' for PCA
## 01:44:00 PCA: 2 components explained 62.53% variance
## 01:44:00 Commencing optimization for 500 epochs, with 179264 positive edges
## 01:44:30 Optimization finished
```

```r
options(ggrepel.max.overlaps = Inf)
metacell_seurat@meta.data$cell_type<-recode(metacell_seurat@meta.data$cell_type, Hepatocyte="Hepatocytes", T_Cells = "T cells", Macrophage="Macrophages", 'NK-like_Cells'="NK-like cells", Erythroid_Cells="Erythroid cells", Mature_B_Cells="Mature B cells", Plasma_Cells="Plasma cells", Portal_endothelial_Cells="Portal endothelial cells")
metacell_plot<-DimPlot(metacell_seurat, reduction = "umap", label = TRUE, repel = TRUE, label.size = 10, group.by = "cell_type", combine = TRUE) + 
  labs(title=paste0("UMAP plot of healthy liver single-cell dataset after metacell transformation", "\n", "(MacParland et al.)"))  +
  theme_classic(base_size = 18) +
  theme(plot.title = element_text(color="black", size=20, face="bold", hjust = 0.5), legend.text=element_text(size=20)) + 
  guides(fill="none", colour = guide_legend(override.aes = list(size=10))) + 
  NoLegend()
plot(metacell_plot)
```

![](markdown_results/unnamed-chunk-4-2.png)<!-- -->

```r
dir.create("./scrna/data/input/")
save(raw_sc_plot, metacell_plot, file="./scrna/data/input/sc_plots.RData")
save(metacell_seurat, file="./scrna/data/input/sc_data.RData")
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
##  [1] BiocParallel_1.28.3 tidyr_1.2.1         tibble_3.1.8       
##  [4] DT_0.24             sp_1.5-0            SeuratObject_4.1.1 
##  [7] Seurat_4.1.0        dplyr_1.0.10        ggplot2_3.3.6      
## [10] cowplot_1.1.1      
## 
## loaded via a namespace (and not attached):
##   [1] Rtsne_0.16            colorspace_2.0-3      deldir_1.0-6         
##   [4] ellipsis_0.3.2        ggridges_0.5.3        spatstat.data_2.2-0  
##   [7] farver_2.1.1          leiden_0.4.2          listenv_0.8.0        
##  [10] ggrepel_0.9.1         RSpectra_0.16-1       fansi_1.0.3          
##  [13] codetools_0.2-18      splines_4.1.2         cachem_1.0.6         
##  [16] knitr_1.40            polyclip_1.10-0       jsonlite_1.8.0       
##  [19] ica_1.0-3             cluster_2.1.3         png_0.1-7            
##  [22] rgeos_0.5-9           uwot_0.1.11           shiny_1.7.2          
##  [25] sctransform_0.3.4     spatstat.sparse_2.1-1 compiler_4.1.2       
##  [28] httr_1.4.4            assertthat_0.2.1      Matrix_1.4-1         
##  [31] fastmap_1.1.0         lazyeval_0.2.2        cli_3.4.0            
##  [34] later_1.2.0           htmltools_0.5.3       tools_4.1.2          
##  [37] igraph_1.3.0          gtable_0.3.1          glue_1.6.2           
##  [40] RANN_2.6.1            reshape2_1.4.4        Rcpp_1.0.9           
##  [43] scattermore_0.8       jquerylib_0.1.4       vctrs_0.4.1          
##  [46] nlme_3.1-159          progressr_0.11.0      lmtest_0.9-40        
##  [49] spatstat.random_2.2-0 xfun_0.32             stringr_1.4.1        
##  [52] globals_0.16.1        mime_0.12             miniUI_0.1.1.1       
##  [55] lifecycle_1.0.1       irlba_2.3.5           goftest_1.2-3        
##  [58] future_1.28.0         MASS_7.3-58.1         zoo_1.8-10           
##  [61] scales_1.2.1          spatstat.core_2.4-4   promises_1.2.0.1     
##  [64] spatstat.utils_2.3-1  parallel_4.1.2        RColorBrewer_1.1-3   
##  [67] yaml_2.3.5            reticulate_1.26       pbapply_1.5-0        
##  [70] gridExtra_2.3         sass_0.4.2            rpart_4.1.16         
##  [73] stringi_1.7.6         highr_0.9             rlang_1.0.5          
##  [76] pkgconfig_2.0.3       matrixStats_0.62.0    evaluate_0.16        
##  [79] lattice_0.20-45       ROCR_1.0-11           purrr_0.3.4          
##  [82] tensor_1.5            labeling_0.4.2        patchwork_1.1.2      
##  [85] htmlwidgets_1.5.4     tidyselect_1.1.2      parallelly_1.32.1    
##  [88] RcppAnnoy_0.0.19      plyr_1.8.7            magrittr_2.0.3       
##  [91] R6_2.5.1              generics_0.1.3        DBI_1.1.3            
##  [94] mgcv_1.8-40           pillar_1.8.1          withr_2.5.0          
##  [97] fitdistrplus_1.1-8    survival_3.4-0        abind_1.4-5          
## [100] future.apply_1.9.1    crayon_1.5.1          KernSmooth_2.23-20   
## [103] utf8_1.2.2            spatstat.geom_2.4-0   plotly_4.10.0        
## [106] rmarkdown_2.16        grid_4.1.2            data.table_1.14.2    
## [109] FNN_1.1.3.1           digest_0.6.29         xtable_1.8-4         
## [112] httpuv_1.6.6          munsell_0.5.0         viridisLite_0.4.1    
## [115] bslib_0.4.0
```
