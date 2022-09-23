---
title: "Preprocessing of plasma cfRNA samples from healthy donors and hepatocellular carcinoma (hcc) patients"
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
library(DESeq2)
library(stringr)
library(dplyr)
library(biomaRt)
library(XML)
library(reutils)
library(tidyr)
library(DT)
```

#  Download the data (GSE142987)
***


```bash
#download the cfRNA count matrix
mkdir -p ./cfrna/data/raw
wget -q -O - https://ftp.ncbi.nlm.nih.gov/geo/series/GSE142nnn/GSE142987/suppl/GSE142987_sample_count_matrix.txt.gz | gunzip -c > ./cfrna/data/raw/cfrna_countm.txt
#download the metadata
wget -q -nv -O ./cfrna/data/raw/cfrna_metadata_init.csv "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/sra-db-be.cgi?rettype=runinfo&term=SRP239389"
```

# Sample metadata construction 
***


```r
cfrna_metadata_init <- read.csv("./cfrna/data/raw/cfrna_metadata_init.csv")
cfrna_metadata_xml <- efetch(c(cfrna_metadata_init$Run), "sra")
cfrna_metadata_xml_cont<-content(cfrna_metadata_xml)
cfrna_metadata_full<-xmlToDataFrame(nodes = getNodeSet(cfrna_metadata_xml_cont, "//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"))
cfrna_metadata_full$samples<-rep(cfrna_metadata_init$Run, each=6)
cfrna_metadata<-cfrna_metadata_full %>% 
  pivot_wider(names_from = TAG, values_from = VALUE)
names(cfrna_metadata)[3]<-"disease_state"
cfrna_metadata$age<-as.numeric(cfrna_metadata$age)
```

#  Cleanup 
***


```r
#load dataset
cfrna_countm <- read.delim("./cfrna/data/raw/cfrna_countm.txt", comment.char="#")
#sort the sample names
cfrna_countm<-cfrna_countm[, str_sort(names(cfrna_countm), numeric = TRUE)]
#transfer gene IDs from column to rownames
rownames(cfrna_countm) <- cfrna_countm$Sample_name
cfrna_countm<-cfrna_countm %>% 
  dplyr::select(-Sample_name)
#shorten and sort new ample names
names(cfrna_countm)<-c(seq(575, 604), seq(540, 574))
cfrna_countm<-cfrna_countm[ , order(names(cfrna_countm))]
#initial dimensions of the data
dim(cfrna_countm)
```

```
## [1] 84656    65
```

```r
#synchronize metadata sample names with count matrix sample names
cfrna_metadata$samples<-str_remove(cfrna_metadata$samples, "SRR10822")
```

#  Normalization and filtering in DESeq2 
***


```r
#input gene expression and metadata into DESeq2 format with the experimental design set to healthy/disease
dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(cfrna_countm),
                                         colData = as.data.frame(cfrna_metadata),
                                         design = ~ disease_state)
#estimate size (normalization) factors
dds <- estimateSizeFactors(dds)
#from the normalized data filter out genes with low expression and keep genes with expression of 5 and higher in at least 90% of samples
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= ncol(cfrna_countm) * 90/100
dds <- dds[idx,]
#get the variance stabilized data
dds<-estimateDispersions(dds)
cfrna_vst<-getVarianceStabilizedData(dds)
#check dimensions 
dim(cfrna_vst)
```

```
## [1] 9200   65
```

```r
#view sample metadata
datatable(cfrna_metadata, caption = 'cfRNA sample characteristics')
```

```{=html}
<div id="htmlwidget-cf921c6df91e7b660e90" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-cf921c6df91e7b660e90">{"x":{"filter":"none","vertical":false,"caption":"<caption>cfRNA sample characteristics<\/caption>","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65"],["540","541","542","543","544","545","546","547","548","549","550","551","552","553","554","555","556","557","558","559","560","561","562","563","564","565","566","567","568","569","570","571","572","573","574","575","576","577","578","579","580","581","582","583","584","585","586","587","588","589","590","591","592","593","594","595","596","597","598","599","600","601","602","603","604"],["liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","liver cancer patient's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma","healthy donor's plasma"],["liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","liver cancer patient","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor","healthy donor"],["male","male","male","female","male","female","male","male","male","male","female","male","male","female","male","male","male","male","male","male","male","male","female","male","male","female","male","male","male","male","male","male","male","male","male","female","female","female","female","female","male","male","female","female","male","female","male","male","male","male","male","female","female","female","male","male","female","male","male","female","male","female","female","female","male"],[43,47,71,57,71,61,57,41,74,71,70,51,43,37,67,71,54,59,43,50,62,74,60,56,47,61,44,80,48,46,52,41,42,51,57,69,65,69,64,58,72,74,66,52,76,62,71,66,57,63,63,54,61,79,54,69,76,63,66,71,65,66,57,69,75],["0","0","0","0","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","B","B","B","B","C","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a","n/a"],["plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma","plasma"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>samples<\/th>\n      <th>source_name<\/th>\n      <th>disease_state<\/th>\n      <th>gender<\/th>\n      <th>age<\/th>\n      <th>Stage<\/th>\n      <th>tissue<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":5},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

#  Getting gene names  
***


```r
httr::set_config(httr::config(ssl_verifypeer = FALSE))
genes<-rownames(cfrna_vst)
#remove version number
genes_clean<-gsub("\\..*","",genes)
#define the dataset to use for conversion
mart <- useEnsembl(dataset="hsapiens_gene_ensembl", biomart="ensembl", version=106) 
#get the hgnc symbols of genes
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("hgnc_symbol", "ensembl_gene_id"), values = genes_clean, mart= mart)
#convert the matrix to dataframe for easier data wrangling
cfrna_vst_df<-as.data.frame(cfrna_vst)
#get the gene IDs into a column
cfrna_vst_df$genes<-genes_clean
#transfer the gene names into the main dataset
cfrna_vst_merge<-merge(cfrna_vst_df, gene_IDs, by.x="genes", by.y="ensembl_gene_id")
#some genes to not have names - substitute the empty spaces with NA
cfrna_vst_merge <- cfrna_vst_merge %>% 
  mutate_all(na_if,"")
#create a new column where we will have the gene names and when names = NA instead we will have gene IDs
cfrna_vst_merge$gene_new<-ifelse(is.na(cfrna_vst_merge$hgnc_symbol), cfrna_vst_merge$genes, cfrna_vst_merge$hgnc_symbol)
#update the rownames
cfrna_vst_merge<-cfrna_vst_merge[!duplicated(cfrna_vst_merge$gene_new), ]
rownames(cfrna_vst_merge) <- cfrna_vst_merge$gene_new
#drop the now unnecessary columns
cfrna_vst_merge <- cfrna_vst_merge %>% 
  dplyr::select(-c(genes, hgnc_symbol, gene_new))
#transform the dataset into a matrix suitable for WGCNA
cfrna_wgcna<- t(as.matrix(cfrna_vst_merge))
#final dimensions of the input data
dim(cfrna_wgcna)
```

```
## [1]   65 9037
```

```r
dir.create("./cfrna/data/input/")
save(cfrna_wgcna, cfrna_metadata, file="./cfrna/data/input/cfrna_input.RData")
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
##  [1] DT_0.24                     tidyr_1.2.1                
##  [3] reutils_0.2.3               XML_3.99-0.9               
##  [5] biomaRt_2.50.0              dplyr_1.0.10               
##  [7] stringr_1.4.1               DESeq2_1.34.0              
##  [9] SummarizedExperiment_1.24.0 Biobase_2.54.0             
## [11] MatrixGenerics_1.6.0        matrixStats_0.62.0         
## [13] GenomicRanges_1.46.1        GenomeInfoDb_1.30.0        
## [15] IRanges_2.28.0              S4Vectors_0.32.3           
## [17] BiocGenerics_0.40.0        
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-7           bit64_4.0.5            filelock_1.0.2        
##  [4] RColorBrewer_1.1-3     progress_1.2.2         httr_1.4.4            
##  [7] tools_4.1.2            bslib_0.4.0            utf8_1.2.2            
## [10] R6_2.5.1               DBI_1.1.3              colorspace_2.0-3      
## [13] tidyselect_1.1.2       prettyunits_1.1.1      bit_4.0.4             
## [16] curl_4.3.2             compiler_4.1.2         cli_3.4.0             
## [19] xml2_1.3.3             DelayedArray_0.20.0    sass_0.4.2            
## [22] scales_1.2.1           genefilter_1.76.0      rappdirs_0.3.3        
## [25] digest_0.6.29          rmarkdown_2.16         XVector_0.34.0        
## [28] pkgconfig_2.0.3        htmltools_0.5.3        dbplyr_2.2.1          
## [31] fastmap_1.1.0          htmlwidgets_1.5.4      rlang_1.0.5           
## [34] RSQLite_2.2.8          jquerylib_0.1.4        generics_0.1.3        
## [37] jsonlite_1.8.0         crosstalk_1.2.0        BiocParallel_1.28.3   
## [40] RCurl_1.98-1.8         magrittr_2.0.3         GenomeInfoDbData_1.2.7
## [43] Matrix_1.4-1           Rcpp_1.0.9             munsell_0.5.0         
## [46] fansi_1.0.3            lifecycle_1.0.1        stringi_1.7.6         
## [49] yaml_2.3.5             zlibbioc_1.40.0        BiocFileCache_2.2.0   
## [52] grid_4.1.2             blob_1.2.3             parallel_4.1.2        
## [55] crayon_1.5.1           lattice_0.20-45        Biostrings_2.62.0     
## [58] splines_4.1.2          annotate_1.72.0        hms_1.1.2             
## [61] KEGGREST_1.34.0        locfit_1.5-9.6         knitr_1.40            
## [64] pillar_1.8.1           geneplotter_1.72.0     glue_1.6.2            
## [67] evaluate_0.16          png_0.1-7              vctrs_0.4.1           
## [70] gtable_0.3.1           purrr_0.3.4            assertthat_0.2.1      
## [73] cachem_1.0.6           ggplot2_3.3.6          xfun_0.32             
## [76] xtable_1.8-4           survival_3.4-0         tibble_3.1.8          
## [79] AnnotationDbi_1.56.1   memoise_2.0.1          ellipsis_0.3.2
```
