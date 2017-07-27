
## Datasets

We used previously published datasets from [Iurlaro et al., 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1001-5).`.bam` files obtained following the details in:

- GEO [GSE77447](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77447)
- ArrayExpress [E-GEOD-77447](http://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-77447/)


## Software

- [samtools v1.2](http://samtools.sourceforge.net/)
- [picard v1.138](https://broadinstitute.github.io/picard/)
- [bedtools v2.25.0](http://bedtools.readthedocs.io/en/latest/)


## Processing

Technical replicates were merged using `samtools merge` and resulting files were sorted and indexed with `samtools sort` and `samtools index` respectively.

Alignment metrics were obtained with `picard.jar CollectAlignmentSummaryMetrics` and `picard.jar picard.jar CollectInsertSizeMetrics`

Whitelist regions for the `mm9` mouse reference genome were prepared from blacklist regions available [here](https://sites.google.com/site/anshulkundaje/projects/blacklists) using `complementBed`.

Filtering was performed with `samtools view` and duplicate reads were marked with `picard.jar MarkDuplicates`


## Peak calling

http://10.20.192.25:3000/projects/20150309-gordon-redbs-mesc/wiki/30-04-2015_Peak_calling
