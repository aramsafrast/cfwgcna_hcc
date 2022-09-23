---
title: "Preprocessing of blood exoRNA samples from hepatocellular carcinoma (hcc) patients"
subtitle: "Data of Li et al. https://doi.org/10.1093/nar/gkx891"
author: "Aram Safrastyan"
date: "18 September, 2022"
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

#  Download and map the data (GSE100207)
* resource/time intensive
***


```bash
#download the exoRNA samples
mkdir -p ./exorna/data/raw/fastq
mkdir -p ./exorna/data/raw/map
mkdir -p ./exorna/data/raw/hs_genome
mkdir -p ./exorna/data/raw/hs_genome/star/
#get human genome fasta sequence
wget -q -O - http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz | gunzip -c > ./exorna/data/raw/hs_genome/hg38.fa
#get human genome gtf annotation
wget -q -O - http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz | gunzip -c > ./exorna/data/raw/hs_genome/hg38.gtf
#make STAR genome index
STAR --runThreadN $threads --runMode genomeGenerate --genomeDir ./exorna/data/raw/hs_genome/star/ --genomeFastaFiles ./exorna/data/raw/hs_genome/hg38.fa --sjdbGTFfile ./exorna/data/raw/hs_genome/hg38.gtf
#download raw sequencing files
for f in {516..536} ; do fasterq-dump SRR5712"$f" -O ./exorna/data/raw/fastq -t ./exorna/data/raw/fastq -e $threads ; pigz -p $threads ./exorna/data/raw/fastq/*fastq ; done
#map with STAR
for f in {516..536} ; do STAR --runMode alignReads --sjdbGTFfile ./exorna/data/raw/hs_genome/hg38.gtf --readFilesIn ./exorna/data/raw/fastq/SRR5712"$f"_1.fastq.gz ./exorna/data/raw/fastq/SRR5712"$f"_2.fastq.gz --readFilesCommand "gunzip -c" --outSAMtype BAM SortedByCoordinate --genomeDir  ./exorna/data/raw/hs_genome/star --outFileNamePrefix ./exorna/data/raw/map/"$f" --runThreadN $threads; done
#generate gene count table 
featureCounts -T $threads -p -s 1 -a ./exorna/data/raw/hs_genome/hg38.gtf -o ./exorna/data/raw/exorna_countm.txt ./exorna/data/raw/map/*bam
#download the metadata
wget -q -nv -O ./exorna/data/raw/exorna_metadata_init.csv "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/sra-db-be.cgi?rettype=runinfo&term=SRP109668"
```

```
## Sep 18 21:02:48 ..... started STAR run
## Sep 18 21:02:48 ... starting to generate Genome files
## Sep 18 21:04:41 ..... processing annotations GTF
## Sep 18 21:05:57 ... starting to sort Suffix Array. This may take a long time...
## Sep 18 21:06:55 ... sorting Suffix Array chunks and saving them to disk...
## Sep 18 22:20:36 ... loading chunks from disk, packing SA...
## Sep 18 22:25:16 ... finished generating suffix array
## Sep 18 22:25:16 ... generating Suffix Array index
## Sep 18 22:32:52 ... completed Suffix Array index
## Sep 18 22:32:53 ..... inserting junctions into the genome indices
## Sep 18 22:39:26 ... writing Genome to disk ...
## Sep 18 22:39:41 ... writing Suffix Array to disk ...
## Sep 18 22:40:51 ... writing SAindex to disk
## Sep 18 22:41:02 ..... finished successfully
## spots read      : 15,337,879
## reads read      : 30,675,758
## reads written   : 30,675,758
## spots read      : 16,318,470
## reads read      : 32,636,940
## reads written   : 32,636,940
## spots read      : 23,907,550
## reads read      : 47,815,100
## reads written   : 47,815,100
## spots read      : 22,620,429
## reads read      : 45,240,858
## reads written   : 45,240,858
## spots read      : 15,610,786
## reads read      : 31,221,572
## reads written   : 31,221,572
## spots read      : 16,331,974
## reads read      : 32,663,948
## reads written   : 32,663,948
## spots read      : 17,229,394
## reads read      : 34,458,788
## reads written   : 34,458,788
## spots read      : 22,611,433
## reads read      : 45,222,866
## reads written   : 45,222,866
## spots read      : 14,594,359
## reads read      : 29,188,718
## reads written   : 29,188,718
## spots read      : 14,770,086
## reads read      : 29,540,172
## reads written   : 29,540,172
## spots read      : 11,693,512
## reads read      : 23,387,024
## reads written   : 23,387,024
## spots read      : 15,834,004
## reads read      : 31,668,008
## reads written   : 31,668,008
## spots read      : 13,734,173
## reads read      : 27,468,346
## reads written   : 27,468,346
## spots read      : 15,696,901
## reads read      : 31,393,802
## reads written   : 31,393,802
## spots read      : 14,616,585
## reads read      : 29,233,170
## reads written   : 29,233,170
## spots read      : 15,863,585
## reads read      : 31,727,170
## reads written   : 31,727,170
## spots read      : 15,644,809
## reads read      : 31,289,618
## reads written   : 31,289,618
## spots read      : 13,697,056
## reads read      : 27,394,112
## reads written   : 27,394,112
## spots read      : 15,930,803
## reads read      : 31,861,606
## reads written   : 31,861,606
## spots read      : 17,289,974
## reads read      : 34,579,948
## reads written   : 34,579,948
## spots read      : 16,615,627
## reads read      : 33,231,254
## reads written   : 33,231,254
## Sep 19 00:27:44 ..... started STAR run
## Sep 19 00:27:44 ..... loading genome
## Sep 19 00:29:26 ..... processing annotations GTF
## Sep 19 00:30:18 ..... inserting junctions into the genome indices
## Sep 19 00:33:10 ..... started mapping
## Sep 19 00:46:35 ..... finished mapping
## Sep 19 00:46:39 ..... started sorting BAM
## Sep 19 00:47:45 ..... finished successfully
## Sep 19 00:47:48 ..... started STAR run
## Sep 19 00:47:48 ..... loading genome
## Sep 19 00:48:27 ..... processing annotations GTF
## Sep 19 00:49:20 ..... inserting junctions into the genome indices
## Sep 19 00:52:11 ..... started mapping
## Sep 19 01:04:33 ..... finished mapping
## Sep 19 01:04:37 ..... started sorting BAM
## Sep 19 01:05:55 ..... finished successfully
## Sep 19 01:05:59 ..... started STAR run
## Sep 19 01:05:59 ..... loading genome
## Sep 19 01:06:40 ..... processing annotations GTF
## Sep 19 01:07:32 ..... inserting junctions into the genome indices
## Sep 19 01:10:21 ..... started mapping
## Sep 19 01:32:45 ..... finished mapping
## Sep 19 01:32:50 ..... started sorting BAM
## Sep 19 01:35:05 ..... finished successfully
## Sep 19 01:35:11 ..... started STAR run
## Sep 19 01:35:11 ..... loading genome
## Sep 19 01:35:50 ..... processing annotations GTF
## Sep 19 01:36:43 ..... inserting junctions into the genome indices
## Sep 19 01:39:35 ..... started mapping
## Sep 19 02:02:24 ..... finished mapping
## Sep 19 02:02:29 ..... started sorting BAM
## Sep 19 02:04:33 ..... finished successfully
## Sep 19 02:04:39 ..... started STAR run
## Sep 19 02:04:39 ..... loading genome
## Sep 19 02:05:18 ..... processing annotations GTF
## Sep 19 02:06:10 ..... inserting junctions into the genome indices
## Sep 19 02:09:03 ..... started mapping
## Sep 19 02:25:45 ..... finished mapping
## Sep 19 02:25:50 ..... started sorting BAM
## Sep 19 02:27:09 ..... finished successfully
## Sep 19 02:27:13 ..... started STAR run
## Sep 19 02:27:13 ..... loading genome
## Sep 19 02:27:54 ..... processing annotations GTF
## Sep 19 02:28:45 ..... inserting junctions into the genome indices
## Sep 19 02:31:37 ..... started mapping
## Sep 19 02:47:18 ..... finished mapping
## Sep 19 02:47:22 ..... started sorting BAM
## Sep 19 02:48:34 ..... finished successfully
## Sep 19 02:48:38 ..... started STAR run
## Sep 19 02:48:38 ..... loading genome
## Sep 19 02:49:17 ..... processing annotations GTF
## Sep 19 02:50:08 ..... inserting junctions into the genome indices
## Sep 19 02:52:59 ..... started mapping
## Sep 19 03:02:52 ..... finished mapping
## Sep 19 03:02:57 ..... started sorting BAM
## Sep 19 03:04:08 ..... finished successfully
## Sep 19 03:04:11 ..... started STAR run
## Sep 19 03:04:11 ..... loading genome
## Sep 19 03:04:51 ..... processing annotations GTF
## Sep 19 03:05:42 ..... inserting junctions into the genome indices
## Sep 19 03:08:31 ..... started mapping
## Sep 19 03:21:31 ..... finished mapping
## Sep 19 03:21:35 ..... started sorting BAM
## Sep 19 03:23:00 ..... finished successfully
## Sep 19 03:23:04 ..... started STAR run
## Sep 19 03:23:04 ..... loading genome
## Sep 19 03:23:43 ..... processing annotations GTF
## Sep 19 03:24:35 ..... inserting junctions into the genome indices
## Sep 19 03:27:27 ..... started mapping
## Sep 19 03:35:23 ..... finished mapping
## Sep 19 03:35:27 ..... started sorting BAM
## Sep 19 03:36:40 ..... finished successfully
## Sep 19 03:36:44 ..... started STAR run
## Sep 19 03:36:44 ..... loading genome
## Sep 19 03:37:24 ..... processing annotations GTF
## Sep 19 03:38:17 ..... inserting junctions into the genome indices
## Sep 19 03:41:07 ..... started mapping
## Sep 19 03:49:21 ..... finished mapping
## Sep 19 03:49:26 ..... started sorting BAM
## Sep 19 03:50:31 ..... finished successfully
## Sep 19 03:50:34 ..... started STAR run
## Sep 19 03:50:34 ..... loading genome
## Sep 19 03:51:13 ..... processing annotations GTF
## Sep 19 03:52:04 ..... inserting junctions into the genome indices
## Sep 19 03:54:54 ..... started mapping
## Sep 19 04:04:35 ..... finished mapping
## Sep 19 04:04:39 ..... started sorting BAM
## Sep 19 04:05:33 ..... finished successfully
## Sep 19 04:05:36 ..... started STAR run
## Sep 19 04:05:36 ..... loading genome
## Sep 19 04:06:17 ..... processing annotations GTF
## Sep 19 04:07:08 ..... inserting junctions into the genome indices
## Sep 19 04:10:01 ..... started mapping
## Sep 19 04:24:54 ..... finished mapping
## Sep 19 04:24:59 ..... started sorting BAM
## Sep 19 04:26:22 ..... finished successfully
## Sep 19 04:26:26 ..... started STAR run
## Sep 19 04:26:26 ..... loading genome
## Sep 19 04:27:06 ..... processing annotations GTF
## Sep 19 04:27:56 ..... inserting junctions into the genome indices
## Sep 19 04:30:49 ..... started mapping
## Sep 19 04:46:27 ..... finished mapping
## Sep 19 04:46:32 ..... started sorting BAM
## Sep 19 04:47:41 ..... finished successfully
## Sep 19 04:47:44 ..... started STAR run
## Sep 19 04:47:44 ..... loading genome
## Sep 19 04:48:22 ..... processing annotations GTF
## Sep 19 04:49:13 ..... inserting junctions into the genome indices
## Sep 19 04:52:03 ..... started mapping
## Sep 19 05:07:12 ..... finished mapping
## Sep 19 05:07:16 ..... started sorting BAM
## Sep 19 05:08:28 ..... finished successfully
## Sep 19 05:08:32 ..... started STAR run
## Sep 19 05:08:32 ..... loading genome
## Sep 19 05:09:10 ..... processing annotations GTF
## Sep 19 05:10:01 ..... inserting junctions into the genome indices
## Sep 19 05:12:52 ..... started mapping
## Sep 19 05:27:23 ..... finished mapping
## Sep 19 05:27:28 ..... started sorting BAM
## Sep 19 05:28:40 ..... finished successfully
## Sep 19 05:28:43 ..... started STAR run
## Sep 19 05:28:43 ..... loading genome
## Sep 19 05:29:22 ..... processing annotations GTF
## Sep 19 05:30:13 ..... inserting junctions into the genome indices
## Sep 19 05:33:05 ..... started mapping
## Sep 19 05:50:08 ..... finished mapping
## Sep 19 05:50:12 ..... started sorting BAM
## Sep 19 05:51:17 ..... finished successfully
## Sep 19 05:51:20 ..... started STAR run
## Sep 19 05:51:20 ..... loading genome
## Sep 19 05:51:58 ..... processing annotations GTF
## Sep 19 05:52:48 ..... inserting junctions into the genome indices
## Sep 19 05:55:38 ..... started mapping
## Sep 19 06:13:27 ..... finished mapping
## Sep 19 06:13:32 ..... started sorting BAM
## Sep 19 06:14:41 ..... finished successfully
## Sep 19 06:14:44 ..... started STAR run
## Sep 19 06:14:44 ..... loading genome
## Sep 19 06:15:22 ..... processing annotations GTF
## Sep 19 06:16:13 ..... inserting junctions into the genome indices
## Sep 19 06:19:03 ..... started mapping
## Sep 19 06:36:52 ..... finished mapping
## Sep 19 06:36:56 ..... started sorting BAM
## Sep 19 06:37:58 ..... finished successfully
## Sep 19 06:38:02 ..... started STAR run
## Sep 19 06:38:02 ..... loading genome
## Sep 19 06:38:40 ..... processing annotations GTF
## Sep 19 06:39:31 ..... inserting junctions into the genome indices
## Sep 19 06:42:23 ..... started mapping
## Sep 19 06:55:44 ..... finished mapping
## Sep 19 06:55:48 ..... started sorting BAM
## Sep 19 06:57:04 ..... finished successfully
## Sep 19 06:57:08 ..... started STAR run
## Sep 19 06:57:08 ..... loading genome
## Sep 19 06:57:47 ..... processing annotations GTF
## Sep 19 06:58:38 ..... inserting junctions into the genome indices
## Sep 19 07:01:30 ..... started mapping
## Sep 19 07:15:04 ..... finished mapping
## Sep 19 07:15:08 ..... started sorting BAM
## Sep 19 07:16:31 ..... finished successfully
## Sep 19 07:16:35 ..... started STAR run
## Sep 19 07:16:35 ..... loading genome
## Sep 19 07:17:16 ..... processing annotations GTF
## Sep 19 07:18:08 ..... inserting junctions into the genome indices
## Sep 19 07:20:59 ..... started mapping
## Sep 19 07:35:07 ..... finished mapping
## Sep 19 07:35:11 ..... started sorting BAM
## Sep 19 07:36:35 ..... finished successfully
## 
##         ==========     _____ _    _ ____  _____  ______          _____  
##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
## 	  v2.0.1
## 
## //========================== featureCounts setting ===========================\\
## ||                                                                            ||
## ||             Input files : 21 BAM files                                     ||
## ||                           o 516Aligned.sortedByCoord.out.bam               ||
## ||                           o 517Aligned.sortedByCoord.out.bam               ||
## ||                           o 518Aligned.sortedByCoord.out.bam               ||
## ||                           o 519Aligned.sortedByCoord.out.bam               ||
## ||                           o 520Aligned.sortedByCoord.out.bam               ||
## ||                           o 521Aligned.sortedByCoord.out.bam               ||
## ||                           o 522Aligned.sortedByCoord.out.bam               ||
## ||                           o 523Aligned.sortedByCoord.out.bam               ||
## ||                           o 524Aligned.sortedByCoord.out.bam               ||
## ||                           o 525Aligned.sortedByCoord.out.bam               ||
## ||                           o 526Aligned.sortedByCoord.out.bam               ||
## ||                           o 527Aligned.sortedByCoord.out.bam               ||
## ||                           o 528Aligned.sortedByCoord.out.bam               ||
## ||                           o 529Aligned.sortedByCoord.out.bam               ||
## ||                           o 530Aligned.sortedByCoord.out.bam               ||
## ||                           o 531Aligned.sortedByCoord.out.bam               ||
## ||                           o 532Aligned.sortedByCoord.out.bam               ||
## ||                           o 533Aligned.sortedByCoord.out.bam               ||
## ||                           o 534Aligned.sortedByCoord.out.bam               ||
## ||                           o 535Aligned.sortedByCoord.out.bam               ||
## ||                           o 536Aligned.sortedByCoord.out.bam               ||
## ||                                                                            ||
## ||             Output file : exorna_countm.txt                                ||
## ||                 Summary : exorna_countm.txt.summary                        ||
## ||              Annotation : hg38.gtf (GTF)                                   ||
## ||      Dir for temp files : ./exorna/data/raw                                ||
## ||                                                                            ||
## ||                 Threads : 15                                               ||
## ||                   Level : meta-feature level                               ||
## ||              Paired-end : yes                                              ||
## ||      Multimapping reads : not counted                                      ||
## || Multi-overlapping reads : not counted                                      ||
## ||   Min overlapping bases : 1                                                ||
## ||                                                                            ||
## ||          Chimeric reads : counted                                          ||
## ||        Both ends mapped : not required                                     ||
## ||                                                                            ||
## \\============================================================================//
## 
## //================================= Running ==================================\\
## ||                                                                            ||
## || Load annotation file hg38.gtf ...                                          ||
## ||    Features : 1499267                                                      ||
## ||    Meta-features : 60708                                                   ||
## ||    Chromosomes/contigs : 47                                                ||
## ||                                                                            ||
## || Process BAM file 516Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 17502155                                             ||
## ||    Successfully assigned alignments : 7878039 (45.0%)                      ||
## ||    Running time : 0.79 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 517Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 19083994                                             ||
## ||    Successfully assigned alignments : 8871487 (46.5%)                      ||
## ||    Running time : 0.91 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 518Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 30572849                                             ||
## ||    Successfully assigned alignments : 11662854 (38.1%)                     ||
## ||    Running time : 1.89 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 519Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 28989327                                             ||
## ||    Successfully assigned alignments : 12153487 (41.9%)                     ||
## ||    Running time : 1.62 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 520Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 19310786                                             ||
## ||    Successfully assigned alignments : 8392431 (43.5%)                      ||
## ||    Running time : 0.83 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 521Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 17022307                                             ||
## ||    Successfully assigned alignments : 9079323 (53.3%)                      ||
## ||    Running time : 0.65 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 522Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 17397001                                             ||
## ||    Successfully assigned alignments : 10116876 (58.2%)                     ||
## ||    Running time : 0.78 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 523Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 19507360                                             ||
## ||    Successfully assigned alignments : 10534583 (54.0%)                     ||
## ||    Running time : 0.93 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 524Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 16525484                                             ||
## ||    Successfully assigned alignments : 8380757 (50.7%)                      ||
## ||    Running time : 0.80 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 525Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 14824896                                             ||
## ||    Successfully assigned alignments : 9053201 (61.1%)                      ||
## ||    Running time : 0.64 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 526Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 12456561                                             ||
## ||    Successfully assigned alignments : 5968878 (47.9%)                      ||
## ||    Running time : 0.54 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 527Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 20326454                                             ||
## ||    Successfully assigned alignments : 6235167 (30.7%)                      ||
## ||    Running time : 1.25 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 528Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 16476129                                             ||
## ||    Successfully assigned alignments : 6913778 (42.0%)                      ||
## ||    Running time : 0.70 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 529Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 16899951                                             ||
## ||    Successfully assigned alignments : 7838025 (46.4%)                      ||
## ||    Running time : 0.72 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 530Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 16594280                                             ||
## ||    Successfully assigned alignments : 8779400 (52.9%)                      ||
## ||    Running time : 0.55 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 531Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 15729485                                             ||
## ||    Successfully assigned alignments : 7917631 (50.3%)                      ||
## ||    Running time : 0.64 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 532Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 15665775                                             ||
## ||    Successfully assigned alignments : 8029551 (51.3%)                      ||
## ||    Running time : 0.71 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 533Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 14232511                                             ||
## ||    Successfully assigned alignments : 7279983 (51.2%)                      ||
## ||    Running time : 0.57 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 534Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 16617075                                             ||
## ||    Successfully assigned alignments : 7894854 (47.5%)                      ||
## ||    Running time : 0.74 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 535Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 19276786                                             ||
## ||    Successfully assigned alignments : 10197671 (52.9%)                     ||
## ||    Running time : 0.87 minutes                                             ||
## ||                                                                            ||
## || Process BAM file 536Aligned.sortedByCoord.out.bam...                       ||
## ||    Strand specific : stranded                                              ||
## ||    Paired-end reads are included.                                          ||
## ||    Total alignments : 18628970                                             ||
## ||    Successfully assigned alignments : 9735088 (52.3%)                      ||
## ||    Running time : 0.69 minutes                                             ||
## ||                                                                            ||
## || Write the final count table.                                               ||
## || Write the read assignment summary.                                         ||
## ||                                                                            ||
## || Summary of counting results can be found in file "./exorna/data/raw/exorn  ||
## || a_countm.txt.summary"                                                      ||
## ||                                                                            ||
## \\============================================================================//
```

# Sample metadata construction 
***


```r
exorna_metadata_init <- read.csv("./exorna/data/raw/exorna_metadata_init.csv")
exorna_metadata_xml <- efetch(c(exorna_metadata_init$Run), "sra")
exorna_metadata_xml_cont<-content(exorna_metadata_xml)
exorna_metadata_full<-xmlToDataFrame(nodes = getNodeSet(exorna_metadata_xml_cont, "//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"))
exorna_metadata_full$samples<-rep(exorna_metadata_init$Run, each=4)
exorna_metadata<-exorna_metadata_full %>% 
  pivot_wider(names_from = TAG, values_from = VALUE)
names(exorna_metadata)[3]<-"disease_state"
```

#  Cleanup 
***


```r
#load dataset
exorna_countm <- read.delim("./exorna/data/raw/exorna_countm.txt", comment.char="#")
exorna_countm <- exorna_countm[, -c(2:6)]
#transfer gene IDs from column to rownames
rownames(exorna_countm) <- exorna_countm$Geneid
exorna_countm<-exorna_countm %>% dplyr::select(-Geneid)
#shorten sample names
names(exorna_countm)<-c(seq(16, 36))
#initial dimensions of the data
dim(exorna_countm)
```

```
## [1] 60708    21
```

```r
#synchronize metadata sample names with count matrix sample names
exorna_metadata$samples<-str_remove(exorna_metadata$samples, "SRR57125")
```

#  Normalization and filtering in DESeq2 
***


```r
#input gene expression and metadata into DESeq2 format with the experimental design set to healthy/disease
dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(exorna_countm),
                                         colData = as.data.frame(exorna_metadata),
                                         design = ~ 1)
#estimate size (normalization) factors
dds <- estimateSizeFactors(dds)
#from the normalized data filter out genes with low expression and keep genes with expression of 5 and higher in at least 90% of samples
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= ncol(exorna_countm) * 90/100
dds <- dds[idx,]
#get the variance stabilized data
dds<-estimateDispersions(dds)
exorna_vst<-getVarianceStabilizedData(dds)
#check dimensions 
dim(exorna_vst)
```

```
## [1] 11485    21
```

```r
#view sample metadata
datatable(exorna_metadata, caption = 'exoRNA sample characteristics')
```

```{=html}
<div id="htmlwidget-98e49d5ef6e90babe4ef" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-98e49d5ef6e90babe4ef">{"x":{"filter":"none","vertical":false,"caption":"<caption>exoRNA sample characteristics<\/caption>","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21"],["16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"],["Blood exosome","Blood exosome","Blood exosome","Blood exosome","Blood exosome","Blood exosome","Blood exosome","Blood exosome","Blood exosome","Blood exosome","Blood exosome","Blood exosome","Blood exosome","Blood exosome","Blood exosome","Blood exosome","Blood exosome","Blood exosome","Blood exosome","Blood exosome","Blood exosome"],["hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma","hepatocellular carcinoma"],["Blood","Blood","Blood","Blood","Blood","Blood","Blood","Blood","Blood","Blood","Blood","Blood","Blood","Blood","Blood","Blood","Blood","Blood","Blood","Blood","Blood"],["exosome","exosome","exosome","exosome","exosome","exosome","exosome","exosome","exosome","exosome","exosome","exosome","exosome","exosome","exosome","exosome","exosome","exosome","exosome","exosome","exosome"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>samples<\/th>\n      <th>source_name<\/th>\n      <th>disease_state<\/th>\n      <th>tissue<\/th>\n      <th>tissue compartment<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

#  Getting gene names  
***


```r
httr::set_config(httr::config(ssl_verifypeer = FALSE))
genes<-rownames(exorna_vst)
#remove version number
genes_clean<-gsub("\\..*","",genes)
#define the dataset to use for conversion
mart <- useEnsembl(dataset="hsapiens_gene_ensembl", biomart="ensembl", version=105) 
#get the hgnc symbols of genes
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("hgnc_symbol", "ensembl_gene_id"), values = genes_clean, mart= mart)
#convert the matrix to dataframe for easier data wrangling
exorna_vst_df<-as.data.frame(exorna_vst)
#get the gene IDs into a column
exorna_vst_df$genes<-genes_clean
#transfer the gene names into the main dataset
exorna_vst_merge<-merge(exorna_vst_df, gene_IDs, by.x="genes", by.y="ensembl_gene_id")
#some genes to not have names - substitute the empty spaces with NA
exorna_vst_merge <- exorna_vst_merge %>% 
  mutate_all(na_if,"")
#create a new column where we will have the gene names and when names = NA instead we will have gene IDs
exorna_vst_merge$gene_new<-ifelse(is.na(exorna_vst_merge$hgnc_symbol), exorna_vst_merge$genes, exorna_vst_merge$hgnc_symbol)
#update the rownames
exorna_vst_merge<-exorna_vst_merge[!duplicated(exorna_vst_merge$gene_new), ]
rownames(exorna_vst_merge) <- exorna_vst_merge$gene_new
#drop the now unnecessary columns
exorna_vst_merge <- exorna_vst_merge %>% 
  dplyr::select(-c(genes, hgnc_symbol, gene_new))
#transform the dataset into a matrix suitable for WGCNA
exorna_wgcna<- t(as.matrix(exorna_vst_merge))
#final dimensions of the input data
dim(exorna_wgcna)
```

```
## [1]    21 11481
```

```r
dir.create("./exorna/data/input/", recursive = TRUE)
save(exorna_wgcna, exorna_metadata, file="./exorna/data/input/exorna_input.RData")
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
