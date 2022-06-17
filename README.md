# WGCNA analysis of cfRNA data from HCC patients

This repistory is aimed at reproducing the results of the paper "Network analysis of hepatocellular carcinoma liquid biopsies augmented by single-cell sequencing data" and is written in R. 

The .Rmd files in the folder *markdown* need to be run in the following order:
1. cfrna_input_clean.Rmd
2. cf_wgcna.Rmd
3. exorna_input_clean.Rmd
4. exo_wgcna.Rmd
5. mp_sc.Rmd (optional as it is resource intensive)

Afterwards, the .Rmd files to generate figures can be run in any desired order. 
