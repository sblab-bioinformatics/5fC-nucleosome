
## Datasets

We used previously published datasets from [Iurlaro et al., 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1001-5). The corresponding `.bam` files were obtained following the details available in:

- GEO [GSE77447](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77447)
- ArrayExpress [E-GEOD-77447](http://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-77447/)


## Software

- [samtools v1.2](http://samtools.sourceforge.net/)
- [picard v1.138](https://broadinstitute.github.io/picard/)
- [bedtools v2.25.0](http://bedtools.readthedocs.io/en/latest/)
- [macs2 v2.1.1.20160226](https://github.com/taoliu/MACS)


## Processing

Technical replicates were merged using `samtools merge` and resulting files were sorted and indexed with `samtools sort` and `samtools index` respectively.

Alignment metrics were obtained with `picard.jar CollectAlignmentSummaryMetrics` and `picard.jar picard.jar CollectInsertSizeMetrics`

Whitelist regions for the `mm9` mouse reference genome were prepared from blacklist regions available [here](https://sites.google.com/site/anshulkundaje/projects/blacklists) using `complementBed`.

Filtering was performed with `samtools view` and duplicate reads were marked with `picard.jar MarkDuplicates`.


## Peak calling

Peaks were called using `macs2 callpeak` with a default q-value threshold of 0.05 and the following parameters:

```bash
macs2 callpeak --keep-dup all -t $input.bam -g 2.5e9 -n $output
```

where `$input.bam` is the input file, `2.5e9` is the effective genome size and `$output` is the output file name. This resulted in `narrowPeak` bed files.

Consensus peaks were obtained using `mergeBed` where a minimum of two peaks overlap by at least 1bp:

```bash
cat *.narrowPeak | sort -k1,1 -k2,2n | awk -v OFS='\t' '{sub("\\..*", "", $4); print $0}' | mergeBed -c 4,4 -o distinct,count_distinct -i - > grm_5fC_merged.bed

awk '$5 > 1' grm_5fC_merged.bed > grm_5fC_consensus.bed
```



