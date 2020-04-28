# Noncoding Project Materials and Code

This repository contains materials and code for reproducing some analyses in noncoding project. For details of project or experimental procedures, please read corresponding manuscript.

## Code

The following code are implemented by Jing Zhang and some of them are (re-)implemented by Shixiang Wang.

- From ICGC mutation data to logistical modeling - [Script 1](pipeline/main.R) - [Script 2](pipeline/2020-p-value-repeat.R) - [Script 3](pipeline/2020-p-value-repeat-02.R)
- From logistical results to final gene list - [Script 1](pipeline/2020-p-value-repeat-03.R) - [Script 2](pipeline/2020-p-value-plotting.R)
- TFBS motif midpoint extraction - [Script](others/TFBS.R)
- Correlation between mutation rate and (epi)genetic features in 1Mb region - [Script 1](pipeline/new_analysis.R) - [Script 2](pipeline/new_analysis2.R)
- Expression correlation between CDC20 and TFBS - [Script](gene_expr_and_association.R)
- ELK4 Chip-Seq signals around CDC promoter region in HELA/HEK293 cell lines - [Script](others/plot_coverage.R) 

Related old repos: [CDC20](https://github.com/XSLiuLab/CDC20), [noncoding](https://github.com/zhangjing1589/noncoding).

## Citation

If you use data, results or conclusion from this work, please cite:

***Pan-cancer noncoding genomic analysis identifies potentially cancer driving CDC20 promoter mutation hotspots, submitted.***


## Abstract

Cancer are primarily caused by somatic mutations in the genomic DNA. Most known cancer driving mutations are located in the protein coding DNA sequence (CDS), and the current exceptions are TERT promoter mutations which are reported to happen in more than 70% of melanoma. Since CDS only occupy less than 2% of the genome sequence, it still remains a mystery if non-coding sequences, which occupy 98% of the human genome, also contain other important cancer driving genetic alterations. Here we performed pan-cancer whole genome mutation analysis to screen for potentially cancer driving mutations in noncoding regions. Recurrent mutations were identified in the promoter of CDC20 gene, which encodes a known key player in cell cycle regulation and cancer progression. In vitro experiments in cell lines support an oncogenic driving function of CDC20 promoter mutations in melanoma progression. In together, our study not only identifies a novel mechanism for CDC20 gene deregulation in human cancers, but also finds novel potentially cancer driving noncoding genetic alterations, with implications for further study of noncoding cancer drivers.

## Steps for Calculating Significance of Recurrent Mutation Regions

### 1. Preprocess ICGC Mutations

We downloaded cancer whole-genome sequencing (WGS) data from International Cancer Genome Consortium (ICGC). In total, there were 4,881 donors, 54,880,488 mutation sites and 59,699,855 mutations before data preprocessing. Nine samples with more than 500,000 mutations were excluded to eliminate ultra-mutated samples. We extracted mutation type “single base substitution” (point mutations) for analysis, and several samples without single base substitution have been removed from analysis. Common human SNP variants were removed from the cancer genome mutation datasets based on 1000 Genomes Project. We also removed the immunoglobulin loci region according to the Ensembl (v75) annotation to avoid bias from immune system-coupled somatic hypermutation. The final mutation data was converted to BED format for subsequent analysis. In total 4859 samples with 47,708,263 mutations are included in downstream analysis. 

### 2. Get Covariates (Genetic or Epigenetic Annotation) for Mutations

We used a variety of annotation features to analyze noncoding mutation rates. These features can be roughly divided into genetic features and epigenetic features. The values of genetic features are correlated with the genomic DNA sequence, and are thus consistent in all tumor types. The value of epigenetic features show variations among cancer types with different tissue origins. These annotation features are described as below.

#### Genetic Features

The genetic annotation features were downloaded from UCSC genome browser database and ENCODE database. 

- Sequence context: We used the 3 base pairs nucleotide motifs centered at the mutated site (1-bp left/right flank motifs of the site). Reverse compliment pairs are combined together, in total there are 32 types of sequence contexts. 
- Genome mappability: This feature refer to the uniqueness of DNA sequence in mapping with reference genome. Genome mappability data was downloaded from UCSC Genome Browser. 
- Recombination rate: Recombination in meiosis help to expand genetic diversity. In somatic cells, DNA lesions can be repaired through recombination between homologous chromosomes. Recombination rate data was downloaded from UCSC Genome Browser. 
- Conservation: We used phastCons method to measure genomic DNA conservation by multiple alignments of 100 vertebrate species in UCSC genome browser. PhastCons estimates the probability that each nucleotide belongs to a conserved element, based on the multiple alignment results.
- Replication timing: We used the ENCODE replication timing data downloaded from the UCSC genome browser. The average wavelet-smoothed signals of repli-seq from 14 cell lines: BJ, GM06990, GM12801, GM12812, GM12813, GM12878, HeLa-S3, HepG2, HUVEC, IMR-90, K-562, MCF-7, NHEK and SK-N-SH were used to assess the genome-wide DNA replication timing.
- GC contents: We used GC content raw data in UCSC genome browser to calculate GC content. The file hg19.gc5Base.txt.gz contains the GC content for 5bp windows across whole genome was downloaded with hgGcPercent. We used this file to recalculate the 1 Mb window GC content across whole genome.
- CpG islands: We used UCSC Genome Table tools to download CpG islands data. The selection criteria is “Mammal”, “Human”, “GRCH37/hg19”, “Regulation”, “CpG Islands”. We calculated the proportion of base pairs in CpG islands in each 1 Mb windows. 
- Promoters: We selected RefSeq-defined human coding genes for analysis. Promoter was defined as the region from 2,500 base pair (bp) upstream to 500 bp downstream from the annotated transcript start site. Pseudogenes are known hot spots for artifacts due to their sequence similarity to their parent genes. In order to avoid potential variant calling bias, partially due to mapping difficulty, we removed the promoters and UTR analyses for pseudogenes. 
- Transcription factor binding sites (TFBS): We used ENCODE union TFBS regions processed by FunSeq2 (http://funseq2.gersteinlab.org/data/2.1.0).
- DNA polymerase II: We used data from ChIP-seq experiments performed by the ENCODE project. The average value of uniform peak signals for 4 cell lines, k562, Mcf10aes, Pbde and Raji were used for analysis. We used these files to recalculate the 1 Mb window DNA polymerase II signal values across whole genome.

#### Epigenetic Features

The epigenetic annotation features were downloaded from Roadmap Epigenomics Project. For pan-cancer analysis, we used the data from integrative analysis of 111 reference human epigenomes. Uniformly processed datasets, integrative analysis products are obtained from this analysis. We used core set of open chromatin and histone modifications consolidated data for further analysis, including chromatin accessibility (DNase-seq), 7 types of histone modifications (H3K4me1, H3K4me3, H3K27me3,  H3K36me3, H3K9me3, H3K27ac and H3K9ac).


### 3. Generate Mutation Specific Probability Model Using Logistic Regression

We used logistic regression model to estimate the background mutation probabilities for each tumor type, patient and position. The position mutation status (0 or 1) is the dependent variable and the genetic and epigenetic features for each position are independent variables. We removed CDS region and immunoglobin loci, and select high-mappable regions from whole genome. In logistic model training, we chose all mutation sites as the positive data, totally 21,857,517 mutation sites. For negative data, we randomly select 50,000,000 sites from whole genome and performed the same select pipeline as mutation site, there were totally 20,534,615 control sites remained.

### 4. Get Region of Interest

We computed the mutation frequency of positive data and extracted the mutation with frequency > 3. Next we constructed 11bp regions by expanding 5bp to the mutaiton position from both left and right. The result regions are then merged.

### 5. Get Region Specific Probabilities

We calculated the mutation probability of each region by summarizing mutation probabilities of all positions in a region.

### 6. Get Poisson Binomial Recurrence Probabilities

We calculated the probabilities of mutation recurrence with Poisson bionomial distribution.

### 7. Get Region Annotations

We extracted start and end position of protein coding genes from reference GTF file (ensembl v75) and queried the upstream 5kb region of gene start position.

### 8. Annotate Recurrent Mutation Regions and Plot

We annotated regions to the closest gene and plotted result gene list.

## CDC20 promoter related database analysis

### CDC20 mRNA expression analysis in pan-cancer

We used the Firebrowse database of Broad Institute to compare the mRNA expression difference between tumor and normal tissues in 37 cancer types. The mRNA expression levels is represented as normalized RSEM (log2).

### Survival analysis

TCGA SKCM patients were selected and divided into two groups, CDC20 mRNA high and CDC20 mRNA low, based on CDC20 mRNA expression level. Kaplan-Meier overall survival curves were compared in these two groups. Log-rank test P value was reported. In ICGC MELA-AU project, we selected patients with mutation occurred in CDC20 locus and nearby regions. Then we divided the patients into two groups, one group with mutation occurring in CDC20 promoter mutation clustered region, another group with mutation occurred in CDC20 locus but not promoter mutation clustered region. Kaplan-Meier overall survival curves were compared in these two groups. 

## Acknowledgement

We thank ICGC and TCGA projects for making cancer genomics data available for analysis. Thank members of Liu lab for helpful discussion. This work was supported in part by the Shanghai Pujiang Program (16PJ1407400), The National Natural Science Foundation of China (31771373), and startup funding from ShanghaiTech University.

## LICENSE

Copyright &copy; 2020 Xue-Song Liu, Jing Zhang, Shixiang Wang

Code and documents of this work are made available for non commercial research purposes only under __Apache License v2__ license, more detail please see [description of license](LICENSE). The code/idea currently may not be used for commercial purposes without explicit written permission after contacting Xue-Song Liu <liuxs@shanghaitech.edu.cn>.

***

**Cancer Biology Group @ShanghaiTech**

**Research group led by Xue-Song Liu in ShanghaiTech. University**
