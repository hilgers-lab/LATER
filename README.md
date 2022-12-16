# ProLoR
## Promoter influence estimation by Long Reads
-------


## Installation 

```
install.packages("devtools")
devtools::install_github("cag1343/ProLoR", build = TRUE, build_vignettes = TRUE)
```

The vignette contains some examples and interpretation of the results of the analysis 
```
library(LORD)
browseVignettes("ProLoR")
```

## Usage

ProLoR estimates transcriptional biases in APA using long read sequencing data 

### Database creation 

First, a database of 5'-3' isoforms is created based on the reference annotation provided. Combinations are computed based on isoform sets TSS and PA sites are merged in a window. This outputs a dataframe with the classification of genes by their TSS and PA site status 


```
linksDatabase <- prepareLinkDatabase(refAnnotation, tss.window=50, tes.window=150)
```

### Counting links 

To account for accurate quantification we develop a counter for long read sequencing data. Aligned reads to the genome are trimmed to their most 5' and 3' end keeping the read identity only reads mapping to both TSS and PA site in the reference, are considered for the analisys. Reads are then summarized in counts per million for further processing. 

```
rtracklayer::import.bed("./genomeAlignments.bed")
countData <- countLinks(alignments, linksDatabase)
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

## Contact

Developer Carlos Alfonso-Gonzalez. For questions or feedback you can contact:

alfonso@ie-freiburg.mpg.de





