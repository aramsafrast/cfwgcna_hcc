---
title: "Generation of Supp. Figs. 1 & 2"
subtitle: "Data of Zhu et al. https://doi.org/10.7150%2Fthno.48206 & Li et al. https://doi.org/10.1093/nar/gkx891"
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
library(dplyr)
library(cowplot)
library(ggplot2)
library(forcats)
library(scales)
library(magick)
library(ragg)
```

# Generation of Supplementary Figure 1
***


```r
load("./cfrna/results/wgcna_main/cfwgcna.RData")
load("./exorna/results/wgcna_main/exowgcna.RData")
supp_fig2<-ggdraw() + draw_plot(cfrna_bar, x=0, y=0, width = .5, height=1) + draw_plot(exorna_bar, x=0.5, y=0, width = .5, height=1) +  draw_plot_label(label = c("A", "B"), size = 28, x = c(0, 0.5), y = c(1, 1), family = "Arial") 
plot(supp_fig2)
```

<img src="markdown_results/supp2-1.png" style="display: block; margin: auto;" />

```r
#save
ggsave(plot=supp_fig2, file="./figures/supp_figures/supp_fig2.png", units = "mm", device = ragg::agg_png, height=80, width=180, scaling = 0.3, limitsize = FALSE)
```

# Generation of Supplementary Figure 2
***


```r
supp_fig1<-ggdraw() + draw_plot(cfrna_param, x=0, y=0.5, width = 1, height=.5) + draw_plot(exorna_param, x=0, y=0, width = 1, height=.5) +  draw_plot_label(label = c("A", "B"), size = 28, x = c(0, 0), y = c(1, 0.5), family = "Arial") 
plot(supp_fig1)
```

<img src="markdown_results/supp-1.png" style="display: block; margin: auto;" />

```r
#save
ggsave(plot=supp_fig1, file="./figures/supp_figures/supp_fig1.png", units = "mm", device = ragg::agg_png, height=150, width=180, scaling = 0.3, limitsize = FALSE)
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
## [1] ragg_1.2.2    magick_2.7.3  scales_1.2.1  forcats_0.5.2 ggplot2_3.3.6
## [6] cowplot_1.1.1 dplyr_1.0.10 
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.9        highr_0.9         pillar_1.8.1      bslib_0.4.0      
##  [5] compiler_4.1.2    jquerylib_0.1.4   tools_4.1.2       digest_0.6.29    
##  [9] jsonlite_1.8.0    evaluate_0.16     lifecycle_1.0.1   tibble_3.1.8     
## [13] gtable_0.3.1      pkgconfig_2.0.3   rlang_1.0.5       cli_3.4.0        
## [17] DBI_1.1.3         yaml_2.3.5        xfun_0.32         fastmap_1.1.0    
## [21] withr_2.5.0       stringr_1.4.1     knitr_1.40        systemfonts_1.0.4
## [25] generics_0.1.3    vctrs_0.4.1       sass_0.4.2        grid_4.1.2       
## [29] tidyselect_1.1.2  glue_1.6.2        R6_2.5.1          textshaping_0.3.6
## [33] fansi_1.0.3       rmarkdown_2.16    farver_2.1.1      purrr_0.3.4      
## [37] magrittr_2.0.3    htmltools_0.5.3   assertthat_0.2.1  colorspace_2.0-3 
## [41] labeling_0.4.2    utf8_1.2.2        stringi_1.7.6     munsell_0.5.0    
## [45] cachem_1.0.6
```
