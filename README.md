# WGCNA analysis of cfRNA data from HCC patients

_____________________________________________

<img src="dag.pdf" align="center" />

_____________________________________________

This repistory is aimed at reproducing the results of the paper [***"Network analysis of hepatocellular carcinoma liquid biopsies augmented by single-cell sequencing data"***](https://doi.org/10.3389/fgene.2022.921195). Front. Genet. 13:921195. The code is written in R. 

For full reproducibility, users are advised to first create a conda enviorenment ([How to install Miniconda]https://docs.conda.io/en/latest/miniconda.html) with the command **conda env create --name cfwgcna_hcc --file=cfwgcna_hcc.yml**. This will install the neccassary tools with the appropriate versions. As some tools were used in original study with different versions, here the newer versions of said tools will be installed. Afterwards, the newly created conda envioremnmtn can be acvited with the command **conda activate cfwgcna_hcc** and the snakemake pipeline to reproduce the results can be run via **snakemake -s snakemake_pipeline --cores X** where "X" is the number of cores one wishes to use. Finally, the figures created with the software "Cytoscape" are not reproduced here and are made already availble but if the user desires the figures can be recreated using the files in the "/cfrna/results/modules/" directory. 

The pipeline has been tested on a Linux system. Using 15 cores the full implementation of the pipeline takes appoximately 12 hours. 

### Contact
For any questions, please contact <aram.safrastyan@uni-jena.de>
