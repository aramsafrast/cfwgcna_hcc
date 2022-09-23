---
title: "Generation of main Fig. 4"
subtitle: "Data of Zhu et al. https://doi.org/10.7150%2Fthno.48206 & MacParland et al. https://doi.org/10.1038/s41467-018-06318-7"
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
library(dplyr)
library(cowplot)
library(ggplot2)
library(forcats)
library(scales)
library(ReactomePA)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ragg)
library(AnnotationDbi)
library(enrichplot)
library(tidyr)
library(Seurat)
library(ggtext)
```

# Generation of Figure 4C, D
***


```r
#for Linux based systems
options(clusterProfiler.download.method = "wget")
#load cfRNA WGCNA results
load("./cfrna/results/wgcna_main/cfwgcna.RData")
#get gene names and transfer to ENTREZID
module_genes<-data.frame( module = cfmergedcolors, gene = colnames(cfrna_wgcna_input)) 
module_genes_entrez<-AnnotationDbi::select(org.Hs.eg.db,
                              keys = module_genes$gene,
                              keytype = "SYMBOL",
                              columns = c("SYMBOL","ENTREZID")) %>% dplyr::left_join(module_genes,by=c("SYMBOL"="gene"))
#generate enrichment plots for module cf-purple
purple_genes<-module_genes_entrez %>% dplyr::filter(module=="cf-purple") %>% 
  na.omit() %>%
  dplyr::select(ENTREZID) %>% 
  unique()
purple_reactome <- enrichPathway(gene= purple_genes$ENTREZID, pvalueCutoff = 0.05)
purple_reactome_df<-fortify(purple_reactome, showCategory = 10) 
purple_reactome_df$p.adjust<- round(-log10(as.numeric(purple_reactome_df$p.adjust)), digits=1)
purple_reactome_plot_old<-enrichplot::cnetplot(purple_reactome,  color_category = "purple", color_gene = "purple", cex_label_category=1, cex_label_gene=1, showCategory = 5) + theme_void(base_size = 24) +   labs(title = "Reactome pathway enrichment analysis", subtitle ="<span style = 'color: red;'>cf-purple</span>") + theme(text=element_text(family="Arial"), plot.margin = unit(c(0, 0, 0, 0), "cm"), plot.subtitle =  element_markdown(size = 25, face="bold", hjust = 0.5), plot.title = element_text(color="black", size=26, face="bold", hjust = 0.5)) + guides(size = guide_legend(reverse=T))
purple_reactome_plot<-ggplot(purple_reactome_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) +  
  geom_point(aes(size = Count, color = p.adjust))  +  
  scale_size_continuous(range = c(6, 10), breaks = pretty_breaks(n=3)) + 
  scale_colour_gradient(name = "-log10(p.adjust)", low="grey", high="black",  breaks=scales::extended_breaks()) +   
  guides(size = guide_legend(reverse=T)) + 
  ylab(NULL) + 
  labs(title = "Reactome pathway enrichment analysis", subtitle ="<span style = 'color: red;'>cf-purple</span>") +  
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 30))  +   
  theme_bw(base_size = 26, base_family = "Arial", base_rect_size = 1, base_line_size = 1) + 
  theme(plot.title = element_text(color="black", size=28, face="bold", hjust = 0.5),
        plot.subtitle =  element_markdown(size = 27, face="bold", hjust = 0.5))
#generate enrichment plots for module cf-yellow
yellow_genes<-module_genes_entrez %>% 
  dplyr::filter(module=="cf-yellow") %>% 
  na.omit() %>%
  dplyr::select(ENTREZID) %>% 
  unique()
yellow_reactome <- enrichPathway(gene= yellow_genes$ENTREZID, pvalueCutoff = 0.05)
yellow_reactome_df<-fortify(yellow_reactome, showCategory = 10) 
yellow_reactome_df$p.adjust<- round(-log10(as.numeric(yellow_reactome_df$p.adjust)), digits=1)
yellow_reactome_plot_old<-enrichplot::cnetplot(yellow_reactome,  color_category = "yellow", color_gene = "yellow", cex_label_category=1, cex_label_gene=1, showCategory = 5) + theme_void(base_size = 24) + 
  labs(title = "Reactome pathway enrichment analysis", subtitle ="<span style = 'color: red;'>cf-yellow</span>") + theme(text=element_text(family="Arial"), plot.margin = unit(c(0, 0, 0, 0), "cm"), plot.subtitle =  element_markdown(size = 25, face="bold", hjust = 0.5), plot.title = element_text(color="black", size=26, face="bold", hjust = 0.5)) + guides(size = guide_legend(reverse=T))
yellow_reactome_plot<-ggplot(yellow_reactome_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) +  
  geom_point(aes(size = Count, color = p.adjust))  +  
  scale_size_continuous(range = c(6, 10), breaks = pretty_breaks(n=3)) + 
  scale_colour_gradient(name = "-log10(p.adjust)", low="grey", high="black",  breaks=scales::extended_breaks()) +   
  guides(size = guide_legend(reverse=T)) + 
  ylab(NULL) + 
  labs(title = "Reactome pathway enrichment analysis", subtitle ="<span style = 'color: red;'>cf-yellow</span>") +  
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 30))  +   
  theme_bw(base_size = 26, base_family = "Arial", base_rect_size = 1, base_line_size = 1) + 
  theme(plot.title = element_text(color="black", size=28, face="bold", hjust = 0.5),
        plot.subtitle =  element_markdown(size = 27, face="bold", hjust = 0.5))
```

# Generation of Figure 4B
***


```r
#load cell specific presevration staistics 
load("./scrna/results/cf_mp/Portal endothelial cells_statsZ.RData")
port<-statsZ_all$Zsummary.pres
load("./scrna/results/cf_mp/Cholangiocytes_statsZ.RData")
chl<-statsZ_all$Zsummary.pres
load("./scrna/results/cf_mp/Hepatocytes_statsZ.RData")
hepato<-statsZ_all$Zsummary.pres
load("./scrna/results/cf_mp/LSECs_statsZ.RData")
lsec<-statsZ_all$Zsummary.pres
load("./scrna/results/cf_mp/Macrophages_statsZ.RData")
macro<-statsZ_all$Zsummary.pres
load("./scrna/results/cf_mp/Mature B cells_statsZ.RData")
b<-statsZ_all$Zsummary.pres
load("./scrna/results/cf_mp/NK-like cells_statsZ.RData")
nk<-statsZ_all$Zsummary.pres
load("./scrna/results/cf_mp/T cells_statsZ.RData")
t<-statsZ_all$Zsummary.pres
load("./scrna/results/cf_mp/Plasma cells_statsZ.RData")
pl<-statsZ_all$Zsummary.pres
load("./scrna/results/cf_mp/Erythroid cells_statsZ.RData")
er<-statsZ_all$Zsummary.pres
#merge, clean and plot the data
cfmp_df<-data.frame(module=rownames(statsZ_all), Portal=port, Hepatocytes=hepato, LSECs=lsec, Macrophages=macro, b_cells=b, t_cells=t, NK_cells=nk, Erythroid=er, Cholangiocytes=chl, Plasma_cells=pl)
cfmp_df$module<-paste0("cf-", "", cfmp_df$module)
rownames(cfmp_df)<-cfmp_df$module 
cfmp_df<-cfmp_df  %>%  
  filter(!module %in% c("cf-gold", "cf-grey")) %>% 
  dplyr::select(-module)
cfmp_df$modules<-rownames(cfmp_df)
cfmp_long<-cfmp_df %>%  
  pivot_longer(!modules, names_to = "cell", values_to = "value") 
cfmp_long$cell <-recode(cfmp_long$cell, t_cells = "T cells", 'NK_cells'="NK-like cells", Erythroid="Erythroid cells", b_cells="Mature B cells", Plasma_cells="Plasma cells", Portal="Portal endothelial cells")
cf_heatmap<-ggplot(cfmp_long, aes(x = cell, y = modules)) +
  geom_tile(color = "black", aes(fill = value)) + theme_classic(base_size = 25) + 
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(2,17), na.value="grey95", name="Zsum") +
  labs(x = "", y = "", fill="cor") +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right") + 
  ggtitle("cfRNA module preservation", subtitle = "scRNA dataset") + 
  theme(
    text=element_text(family="Arial"),
    legend.position="left",
    legend.title.align=0.5,
    legend.key = element_rect(fill = "lightblue", color = NA),
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1,"cm"), 
    plot.title = element_text(color="black", size=26, face="bold", hjust = 0.5),
    plot.subtitle = element_text(color="black", size=25, face="bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, vjust=1, hjust=1),
    axis.text.y = element_text(angle = 0, hjust = 1, colour = c("black", "black", "black", "black", "red", "black", "black", "red"))) 
```

# Generation of Figure 4
***


```r
#load scRNA UMAP plot
load("./scrna/data/input/sc_plots.RData")
#combine the plots and save
options(ggrepel.max.overlaps = Inf)
metacell_plot<-metacell_plot + 
  labs(title="UMAP plot of liver single-cell dataset", subtitle = "Metacell transformed (MacParland et al.)")  +
  theme_classic(base_size = 24, base_family = "Arial") +
  theme(text=element_text(family="Arial"), plot.title = element_text(color="black", size=26, face="bold", hjust = 0.5), plot.subtitle = element_text(color="black", size=25, face="bold", hjust = 0.5), legend.text=element_text(size=7)) + 
  guides(fill="none", colour = guide_legend(override.aes = list(size=10))) + 
  NoLegend()
fig4<-ggdraw() + 
  draw_plot(metacell_plot, x=0, y=0.5, width=.5, height = .5) + 
  draw_plot(cf_heatmap, x=0.5, y=0.5, width=.5, height = .5) + 
  draw_plot(yellow_reactome_plot, x=0, y=0, height=.5, width = .5) + 
  draw_plot(purple_reactome_plot, x=0.5, y=0, height=0.5, width = .5) + 
  draw_plot_label(label = c("A", "B", "C", "D"), size = 38, x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.5, 0.5), family = "Arial") 
plot(fig4)
```

<img src="markdown_results/fig4-1.png" style="display: block; margin: auto;" />

```r
ggsave(plot=fig4, file="./figures/main_figures/fig4.png", units = "mm", device = ragg::agg_png, height=120, width=180, scaling = 0.3, limitsize = FALSE)
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
##  [1] ggtext_0.1.1          sp_1.5-0              SeuratObject_4.1.1   
##  [4] Seurat_4.1.0          tidyr_1.2.1           enrichplot_1.14.1    
##  [7] ragg_1.2.2            clusterProfiler_4.2.0 org.Hs.eg.db_3.14.0  
## [10] AnnotationDbi_1.56.1  IRanges_2.28.0        S4Vectors_0.32.3     
## [13] Biobase_2.54.0        BiocGenerics_0.40.0   ReactomePA_1.38.0    
## [16] scales_1.2.1          forcats_0.5.2         ggplot2_3.3.6        
## [19] cowplot_1.1.1         dplyr_1.0.10         
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.2.2             reticulate_1.26        tidyselect_1.1.2      
##   [4] RSQLite_2.2.8          htmlwidgets_1.5.4      grid_4.1.2            
##   [7] BiocParallel_1.28.3    Rtsne_0.16             scatterpie_0.1.8      
##  [10] munsell_0.5.0          codetools_0.2-18       ica_1.0-3             
##  [13] future_1.28.0          miniUI_0.1.1.1         withr_2.5.0           
##  [16] spatstat.random_2.2-0  colorspace_2.0-3       GOSemSim_2.20.0       
##  [19] progressr_0.11.0       highr_0.9              knitr_1.40            
##  [22] ROCR_1.0-11            tensor_1.5             DOSE_3.20.0           
##  [25] listenv_0.8.0          labeling_0.4.2         GenomeInfoDbData_1.2.7
##  [28] polyclip_1.10-0        bit64_4.0.5            farver_2.1.1          
##  [31] downloader_0.4         parallelly_1.32.1      vctrs_0.4.1           
##  [34] treeio_1.18.0          generics_0.1.3         xfun_0.32             
##  [37] markdown_1.1           R6_2.5.1               GenomeInfoDb_1.30.0   
##  [40] graphlayouts_0.8.1     bitops_1.0-7           spatstat.utils_2.3-1  
##  [43] cachem_1.0.6           fgsea_1.20.0           gridGraphics_0.5-1    
##  [46] assertthat_0.2.1       promises_1.2.0.1       ggraph_2.0.6          
##  [49] rgeos_0.5-9            gtable_0.3.1           globals_0.16.1        
##  [52] goftest_1.2-3          tidygraph_1.2.2        rlang_1.0.5           
##  [55] systemfonts_1.0.4      splines_4.1.2          lazyeval_0.2.2        
##  [58] spatstat.geom_2.4-0    checkmate_2.1.0        yaml_2.3.5            
##  [61] reshape2_1.4.4         abind_1.4-5            backports_1.4.1       
##  [64] httpuv_1.6.6           qvalue_2.26.0          gridtext_0.1.4        
##  [67] tools_4.1.2            ggplotify_0.1.0        ellipsis_0.3.2        
##  [70] spatstat.core_2.4-4    jquerylib_0.1.4        RColorBrewer_1.1-3    
##  [73] ggridges_0.5.3         Rcpp_1.0.9             plyr_1.8.7            
##  [76] zlibbioc_1.40.0        purrr_0.3.4            RCurl_1.98-1.8        
##  [79] rpart_4.1.16           deldir_1.0-6           pbapply_1.5-0         
##  [82] viridis_0.6.2          zoo_1.8-10             ggrepel_0.9.1         
##  [85] cluster_2.1.3          magrittr_2.0.3         data.table_1.14.2     
##  [88] scattermore_0.8        DO.db_2.9              lmtest_0.9-40         
##  [91] RANN_2.6.1             reactome.db_1.77.0     fitdistrplus_1.1-8    
##  [94] matrixStats_0.62.0     patchwork_1.1.2        mime_0.12             
##  [97] evaluate_0.16          xtable_1.8-4           gridExtra_2.3         
## [100] compiler_4.1.2         tibble_3.1.8           KernSmooth_2.23-20    
## [103] crayon_1.5.1           shadowtext_0.1.2       htmltools_0.5.3       
## [106] mgcv_1.8-40            ggfun_0.0.7            later_1.2.0           
## [109] aplot_0.1.7            DBI_1.1.3              tweenr_2.0.2          
## [112] MASS_7.3-58.1          rappdirs_0.3.3         Matrix_1.4-1          
## [115] cli_3.4.0              parallel_4.1.2         igraph_1.3.0          
## [118] pkgconfig_2.0.3        plotly_4.10.0          spatstat.sparse_2.1-1 
## [121] xml2_1.3.3             ggtree_3.2.0           bslib_0.4.0           
## [124] XVector_0.34.0         yulab.utils_0.0.5      stringr_1.4.1         
## [127] digest_0.6.29          sctransform_0.3.4      RcppAnnoy_0.0.19      
## [130] graph_1.72.0           spatstat.data_2.2-0    Biostrings_2.62.0     
## [133] rmarkdown_2.16         leiden_0.4.2           fastmatch_1.1-3       
## [136] tidytree_0.4.0         uwot_0.1.11            shiny_1.7.2           
## [139] graphite_1.40.0        lifecycle_1.0.1        nlme_3.1-159          
## [142] jsonlite_1.8.0         viridisLite_0.4.1      fansi_1.0.3           
## [145] pillar_1.8.1           lattice_0.20-45        KEGGREST_1.34.0       
## [148] fastmap_1.1.0          httr_1.4.4             survival_3.4-0        
## [151] GO.db_3.14.0           glue_1.6.2             png_0.1-7             
## [154] bit_4.0.4              ggforce_0.3.4          stringi_1.7.6         
## [157] sass_0.4.2             blob_1.2.3             textshaping_0.3.6     
## [160] memoise_2.0.1          irlba_2.3.5            future.apply_1.9.1    
## [163] ape_5.6-2
```
