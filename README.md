<!-- badges: start -->

![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/hilgers-lab/ProLoR)
[![Maintained?](https://img.shields.io/badge/Maintained%3F-Yes-brightgreen)](https://github.com/hilgers-lab/ProLoR/graphs/contributors)
[![Install](https://img.shields.io/badge/Install-Github-brightgreen)](#installation)
[![Downloads](https://img.shields.io/github/downloads/hilgers-lab/ProLoR/total)]()
![GitHub](https://img.shields.io/github/license/hilgers-lab/ProLoR)
<!-- badges: end -->

# ProLoR
## Promoter influence estimation by Long Reads
-------


## Installation 

```
install.packages("devtools")
devtools::install_github("hilgers-lab/ProLoR", build = TRUE, build_vignettes = TRUE)
```

The vignette contains some examples and interpretation of the results of the analysis 
```
library(ProLoR)
vignette("ProLoR")
```

## Usage

ProLoR estimates transcriptional biases in APA using long read sequencing data 

# Input data: 
  * Genome Alignment bam files [minimap2](https://github.com/lh3/minimap2) using parameters `minimap2 -ax splice -u f annotation/genome.fa long_read.fastq.gz | samtools sort -@ 4 -o output.bam - samtools index output.bam`
  * Reference annotation in gtf format. Example file [here](https://github.com/hilgers-lab/ProLoR/blob/master/inst/exdata/dm6.annot.gtf.gz) 

### Database creation 

First, a database of 5'-3' isoforms is created based on the reference annotation provided. Combinations are computed based on isoform sets TSS and PA sites are merged in a window. This outputs a dataframe with the classification of genes by their TSS and PA site status 


```
annot_path <- system.file("exdata/dm6.annot.gtf.gz", package="ProLoR")
refAnnotation <- rtracklayer::import.gff(annot_path)
linksDatabase <- prepareLinkDatabase(refAnnotation, tss.window=50, tes.window=150)
```

### Counting links 

To account for accurate quantification we develop a counter for long read sequencing data. Aligned reads to the genome are trimmed to their most 5' and 3' end keeping the read identity only reads mapping to both TSS and PA site in the reference, are considered for the analisys. Reads are then summarized in counts per million for further processing. 

```
bamPath <- system.file("exdata/testBam.bam", package = 'ProLoR')
countData <- countLinks(bamPath, linksDatabase)
```


### Estimating promoter dominance 

Promoter dominance estimates are calculated as perfomed in (Alfonso-Gonzalez, et al., 2022). This function outputs per promoter biases in expression of a given 3'end of the gene. 

```
promoterContributionEstimates <- calculatePromoterDominance(countData, linksDatabase$pairsDatabase)
```


### Estimating transcriptional biases 

Transcriptional biases are calculated by estimating using the joint frequencies of TSS-PA site combinations per gene. Coupling events per gene are estimated using multinomial testing using chi-square. Statistical testing is also available with `fisher.test()` using method="fisher". 

```
biasGenes <- estimateTranscriptionalBias(countData, linksDatabase$pairsDatabase, method="fisher")
```
# Release 

Initial Release 0.1.0

Release date: 22th Dec 2022
This release corresponds to the ProLoR version used by Alfonso-Gonzalez et al. manuscript

## Contact

Developer Carlos Alfonso-Gonzalez. For questions or feedback you can contact:

alfonso@ie-freiburg.mpg.de





