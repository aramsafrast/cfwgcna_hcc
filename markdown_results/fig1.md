---
knit: (function(inputFile, encoding) {
  rmarkdown::render(inputFile, encoding = encoding, output_dir = "./markdown_results") })
title: "Generation of main Fig. 1"
subtitle: "Data of Zhu et al. https://doi.org/10.7150%2Fthno.48206"
author: "Aram Safrastyan"
date: "24 Juni, 2022"
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
library(WGCNA)
library(dplyr)
library(cowplot)
library(ggplot2)
library(tibble)
library(ggrepel)
library(tidyr)
library(ragg)
library(scales)
library(forcats)
library(stringr)
library(openxlsx)
library(DT)
```

# Generation of Fig. 1A
***


```r
load("./cfrna/results/wgcna_main/cfwgcna.RData")
#get the gene and sample numbers of gene expression input
nGenes <- ncol(cfrna_wgcna_input)
nSamples <- nrow(cfrna_wgcna_input)
#get and order module eigengenes 
MEs0 <- moduleEigengenes(cfrna_wgcna_input, cfmergedcolors)$eigengenes
MEs <- orderMEs(MEs0)
#remove grey module which contains unassigned genes
MEs<- removeGreyME(MEs, greyMEName = paste(moduleColor.getMEprefix(), "grey", sep=""))
#correlate traits with modules 
moduleTraitCor <- cor(MEs, cfrna_meta_input, use = "p")
#get the p values
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
#FDR correction of p values
moduleTraitPvalue[ ,1] <- p.adjust(moduleTraitPvalue[ ,1], method = "fdr")
moduleTraitPvalue[ ,2] <- p.adjust(moduleTraitPvalue[ ,2], method = "fdr")
moduleTraitPvalue[ ,3] <- p.adjust(moduleTraitPvalue[ ,3], method = "fdr")
#to plot the results in ggplot some data wrangling needs to be done
module_trait_df<-as.data.frame(moduleTraitCor)
module_trait_df$modules<-gsub("ME", "", rownames(module_trait_df))
module_trait_long<-module_trait_df %>%  
  pivot_longer(!modules, names_to = "cor", values_to = "value") %>% 
  mutate(cor = fct_recode(cor, "disease state" = "disease_state"))
module_p_df<-as.data.frame(moduleTraitPvalue)
module_p_df$modules<-gsub("ME", "", rownames(module_p_df))
module_p_long<-module_p_df %>%  
  pivot_longer(!modules, names_to = "cor", values_to = "value") %>% 
  mutate(cor = fct_recode(cor, "disease state" = "disease_state"))
module_trait_long$modules<-factor(module_trait_long$modules, levels=c("cf-blue", "cf-purple", "cf-brown", "cf-red", "cf-magenta", "cf-black", "cf-yellow", "cf-turquoise"))
heatmap_core <- ggplot(module_trait_long, aes(x = cor, y = modules)) +
  geom_tile(color = "black", aes(fill = value)) + 
  scale_fill_distiller(palette = "RdBu") +
  labs(x = "", y = "", fill="cor") +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0)) + 
  theme(
    legend.title.align=0.5,
    legend.key = element_rect(fill = "lightblue", color = NA),
    legend.key.size = unit(2.5, "cm"),
    legend.key.width = unit(1.5,"cm") 
  ) 
heatmap_full <- heatmap_core + 
  geom_point(data = module_p_long, aes(x=cor, y=modules, size=-log10(value)), alpha=0.5) + 
  scale_size("-log10(p.adjust)", range=c(5, 20), breaks=c(5, 10, 25)) + 
  guides(size = guide_legend(reverse=T)) +   
  theme_bw(base_size = 17, base_family = "Arial", base_rect_size = 1, base_line_size = 1) +
  theme(plot.title = element_text(color="black", size=18, face="bold", hjust = 0.5), plot.subtitle = element_text(color="black", size=17, face="bold", hjust = 0.5)) +
  theme(legend.key=element_blank(), legend.key.height = unit(1, "cm"), legend.key.width = unit(1, "cm")) + 
  ggtitle("Module-trait relationships for cell-free RNA dataset", subtitle = ("(Zhu et al.)"))
print(heatmap_full)
```

<img src="./markdown_results/corr-1.png" style="display: block; margin: auto;" />

# Generation of Fig. 1B,C
***


```r
#extract the cancer stage information
condition <- as.data.frame(cfrna_meta_input$disease_state)
names(condition) <- "Condition"
#names (colors) of the modules
modNames <- substring(names(MEs), 3)
#calculate the gene module membership and its significance 
geneModuleMembership <- as.data.frame(cor(cfrna_wgcna_input, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")
#calculate the correlation and significance of gene-trait relationship
geneTraitSignificance <- as.data.frame(cor(cfrna_wgcna_input, condition, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(condition), sep="")
names(GSPvalue) <- paste("p.GS.", names(condition), sep="")
geneModuleMembership$trait<- geneTraitSignificance$GS.Condition
#choose the module with the highest negative correlation with disease condition
module="cf-blue"
column <- match(module, modNames)
moduleGenes <- cfmergedcolors==module
#plot the relationship between gene correlation with disease condition and gene module membership
blue<-geneModuleMembership %>%
  filter(moduleGenes) %>%
  rownames_to_column(., var = "genes") %>% {
  ggplot(., aes(x=`MMcf-blue`, y=trait))  +
  geom_point(size=2, color="blue", alpha=0.15)  +
  geom_smooth(method="lm", color="red") + 
  scale_x_continuous(limits = c(-1, 1))  + 
  scale_y_continuous(limits = c(-1, 1))  +  
  labs(x = "Module membership in cf-blue module", y = paste0("Gene-trait significance", "\n", "for disease state"), title="Module membership vs. gene significance", subtitle = "module cf-blue")  + 
       theme_classic(base_size = 17, base_family = "Arial") +
      theme(plot.title = element_text(color="black", size=18, face="bold", hjust = 0.5),
      plot.subtitle = element_text(color="black", size=17, face="bold", hjust = 0.5)) 
  }
blue
```

<img src="./markdown_results/mm-1.png" style="display: block; margin: auto;" />

```r
#choose the module with the highest positive correlation with disease condition
module="cf-turquoise"
column <- match(module, modNames)
moduleGenes <- cfmergedcolors==module
turq<-geneModuleMembership %>%
  filter(moduleGenes) %>%
  rownames_to_column(., var = "genes") %>% {
  ggplot(., aes(x=`MMcf-turquoise`, y=trait))  +
  geom_point(size=2, color="turquoise", alpha=0.15)  +
  geom_smooth(method="lm", color="red") + 
  scale_x_continuous(limits = c(-1, 1))  + 
  scale_y_continuous(limits = c(-1, 1))  + 
  labs(x = "Module membership in cf-turquoise module", y = paste0("Gene-trait significance", "\n", "for disease state"), title="Module membership vs. gene significance",  subtitle = "module cf-turquoise")  +
     theme_classic(base_size = 17, base_family = "Arial") +
      theme(plot.title = element_text(color="black", size=18, face="bold", hjust = 0.5),
       plot.subtitle = element_text(color="black", size=17, face="bold", hjust = 0.5)) 
  }
turq
```

<img src="./markdown_results/mm-2.png" style="display: block; margin: auto;" />

# Generation of Fig. 1
***


```r
fig1<-ggdraw() +
  draw_plot(heatmap_full, x = 0, y = 0, height=1, width=.55) +
  draw_plot(blue, x = .55, y = 0.5, width=.45, height = .5) +
  draw_plot(turq, x = .55, y = 0, width =.45, height = .5) +
  draw_plot_label(label = c("A", "B", "C"), size = 24,
                  x = c(0, 0.55, 0.55), y = c(0.98, 0.98, 0.48), family = "Arial")
fig1
```

<img src="./markdown_results/fig1-1.png" style="display: block; margin: auto;" />

```r
dir.create("./figures/")
dir.create("./figures/main_figures/")
ggsave(plot=fig1, file="./figures/main_figures/fig1.png", units = "mm", device = ragg::agg_png, height=100, width=180, scaling = 0.445, limitsize = FALSE)
```


```r
sessionInfo()
```

```
## R version 4.1.2 (2021-11-01)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 22.04 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0
## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] DT_0.23               openxlsx_4.2.5        stringr_1.4.0        
##  [4] forcats_0.5.1         scales_1.2.0          ragg_1.2.2           
##  [7] tidyr_1.2.0           ggrepel_0.9.1         tibble_3.1.7         
## [10] ggplot2_3.3.6         cowplot_1.1.1         dplyr_1.0.9          
## [13] WGCNA_1.71            fastcluster_1.2.3     dynamicTreeCut_1.63-1
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-155           bitops_1.0-7           matrixStats_0.62.0    
##  [4] bit64_4.0.5            doParallel_1.0.17      RColorBrewer_1.1-3    
##  [7] httr_1.4.3             GenomeInfoDb_1.30.1    tools_4.1.2           
## [10] backports_1.4.1        bslib_0.3.1            utf8_1.2.2            
## [13] R6_2.5.1               rpart_4.1.16           mgcv_1.8-39           
## [16] Hmisc_4.7-0            DBI_1.1.2              BiocGenerics_0.40.0   
## [19] colorspace_2.0-3       nnet_7.3-17            withr_2.5.0           
## [22] tidyselect_1.1.2       gridExtra_2.3          preprocessCore_1.56.0 
## [25] bit_4.0.4              compiler_4.1.2         textshaping_0.3.6     
## [28] cli_3.3.0              Biobase_2.54.0         htmlTable_2.4.0       
## [31] labeling_0.4.2         sass_0.4.1             checkmate_2.1.0       
## [34] systemfonts_1.0.4      digest_0.6.29          foreign_0.8-82        
## [37] rmarkdown_2.14         XVector_0.34.0         base64enc_0.1-3       
## [40] jpeg_0.1-9             pkgconfig_2.0.3        htmltools_0.5.2       
## [43] highr_0.9              fastmap_1.1.0          htmlwidgets_1.5.4     
## [46] rlang_1.0.2            impute_1.68.0          rstudioapi_0.13       
## [49] RSQLite_2.2.14         farver_2.1.0           jquerylib_0.1.4       
## [52] generics_0.1.2         jsonlite_1.8.0         zip_2.2.0             
## [55] RCurl_1.98-1.7         magrittr_2.0.3         GO.db_3.14.0          
## [58] GenomeInfoDbData_1.2.7 Formula_1.2-4          Matrix_1.4-0          
## [61] Rcpp_1.0.8.3           munsell_0.5.0          S4Vectors_0.32.4      
## [64] fansi_1.0.3            lifecycle_1.0.1        stringi_1.7.6         
## [67] yaml_2.3.5             zlibbioc_1.40.0        grid_4.1.2            
## [70] blob_1.2.3             parallel_4.1.2         crayon_1.5.1          
## [73] lattice_0.20-45        Biostrings_2.62.0      splines_4.1.2         
## [76] KEGGREST_1.34.0        knitr_1.39             pillar_1.7.0          
## [79] codetools_0.2-18       stats4_4.1.2           glue_1.6.2            
## [82] evaluate_0.15          latticeExtra_0.6-29    data.table_1.14.2     
## [85] png_0.1-7              vctrs_0.4.1            foreach_1.5.2         
## [88] gtable_0.3.0           purrr_0.3.4            assertthat_0.2.1      
## [91] cachem_1.0.6           xfun_0.31              survival_3.2-13       
## [94] iterators_1.0.14       AnnotationDbi_1.56.2   memoise_2.0.1         
## [97] IRanges_2.28.0         cluster_2.1.2          ellipsis_0.3.2
```
