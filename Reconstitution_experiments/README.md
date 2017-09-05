
## Software

- [fastqc v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [bwa v0.7.12-r1044](http://bio-bwa.sourceforge.net/)
- [samtools v1.2](http://samtools.sourceforge.net/)
- [BamUtil v1.0.13](https://genome.sph.umich.edu/wiki/BamUtil)
- [picard v1.138](https://broadinstitute.github.io/picard/)


## Processing

This is mouse brain tissue (11.5 days of development) where genomic DNA was extracted and nucleosome organisation was reconstituted. There were four libraries: 

- two reduced libraries with all 5fC converted to 5hmC (A012 and A013)
- two normal libraries with normal 5fC (A004 and A006)

Sequencing quality was explored using `fastqc` and looked acceptable for further analysis. Illumina sequencing adaptors were removed from the paired end sequencing data of mouse hindbrain and heart using `cutadapt`:

```bash
cutadapt -m 10 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o R1_trimmed.fastq.gz -p R2_trimmed.fastq.gz R1.fastq.gz R2.fastq.gz > output_statistics.txt
```

where `R1.fastq.gz` and `R2.fastq.gz` are the R1 and R2 paired end raw sequencing files with the `_trimmed` counterparts obtained as output. Alignments to the mouse reference genome `mm9` were generated using `bwa mem`:

```bash
bwa mem -M -t 8 $ref $fq1 $fq2
```

where `$ref` is the fasta reference genome file, and `$fq1` and `$fq2` correspond to the trimmed R1 and R2 sequencing files respectively.

Unmapped, not primary aligned and supplementary reads were filtered from the resulting alignment `.bam` files using `samtools view`. Reads aligning with low mapping quality or aligning to [blacklisted genomic regions](https://sites.google.com/site/anshulkundaje/projects/blacklists) of the `mm9` reference genome were also filtered out.

```bash
samtools view -@8 -S -u -F2820 -q 5 -L mm9-whitelist.bed input.bam | samtools sort -@8 - output.bam
```

where `mm9-whitelist.bed` contains not blacklisted genomic regions and `input.bam` is the alignment file obtained with `bwa mem`. Resulting filtered alignment files were sorted by coordinate using `samtools sort`. Overlapping read pairs were clipped using `BamUtil` tools so they do not overlap:

```bash
bam clipOverlap --in input.bam --out output.bam --stats --storeOrig XC
```

Sequencing duplicates were marked using `picard`:

```bash
java -Xmx3G -jar picard.jar MarkDuplicates I=input.bam O=output.bam M=log.markdup.txt
```

were `input.bam` is the merged and sorted alignment file, with duplicates marked as `output.bam` and processing log details written to `log.markdup.txt`. Resulting files were indexed using `samtools index`.

Details about the insert size were obtained using `picard` too:

```bash
java -Xmx2G -jar picard.jar CollectInsertSizeMetrics I=input.bam O=InsertSize.txt H=InsertSize.pdf AS=true VALIDATION_STRINGENCY=SILENT
```




## Generate sequences from 5fC-sites flanking regions

TODO --> how did we generate `TDG_KO.Hindbrain.cns4.unique.end.fasta` and the rest? 

## Computational prediction of nucleosome occupancy

The code used, along with its usage instructions, is described [here](Comp_prediction_nucleosome_occupancy/)
