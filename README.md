# Post-GWAS

## Enrichment Analysis
Test if significant SNPs are enriched in given genomic features (e.g., highly conserved elements, UTRs, differentially expressed genes) through randomization of genomic positions of the genomic features. To speed up the randomization, we can ignore gaps and merge all chromosomes into continuous one.

To run the script, type the command below. With this command, we would get an output file named out_sd.1-1000.txt. 
```
perl fast_enrichment.pl umd3.1.chrom.bed sig_snps.bed sd.bed out_sd 1 1000
```
All input files should be in BED format, and only the first three columns are required.

When genomic gaps should be kept in enrichment analysis, try to use another script of mine, [sd-analysis/association_analysis](https://github.com/jiang18/sd-analysis/tree/master/association_analysis).

### Citation
Jiang, Jicai, et al. "Global copy number analyses by next generation sequencing provide insight into pig genome variation." Bmc Genomics 15.1 (2014): 593.
