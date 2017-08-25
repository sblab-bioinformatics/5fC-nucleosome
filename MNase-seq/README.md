
## Software

- [fastqc v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [bwa v0.7.12-r1044](http://bio-bwa.sourceforge.net/)
- [samtools v1.2](http://samtools.sourceforge.net/)
- [picard v1.138](https://broadinstitute.github.io/picard/)
- [bedtools v2.25.0](http://bedtools.readthedocs.io/en/latest/)


## Processing

Sequencing quality was explored using `fastqc` and looked acceptable for further analysis. Illumina sequencing adaptors were removed from the paired end sequencing data of mouse hindbrain and heart using `cutadapt`:

```bash
cutadapt -m 10 -O 3 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o R1_trimmed.fastq.gz -p R2_trimmed.fastq.gz R1.fastq.gz R2.fastq.gz
```

where `R1.fastq.gz` and `R2.fastq.gz` are the R1 and R2 paired end raw sequencing files with the `_trimmed` counterparts obtained as output. Alignments to the mouse reference genome `mm9` were generated using `bwa mem`:

```bash
bwa mem -M -t 8 $ref $fq1 $fq2
```

where `$ref` is the fasta reference genome file, and `$fq1` and `$fq2` correspond to the trimmed R1 and R2 sequencing files respectively.

Unmapped, not primary aligned and supplementary reads were filtered from the resulting alignment `.bam` files using `samtools view`. Reads aligning with low mapping quality or aligning to [blacklisted genomic regions](https://sites.google.com/site/anshulkundaje/projects/blacklists) of the `mm9` reference genome were also filtered out.

```bash
samtools view -@8 -S -u -F2820 -q 5 -L mm9-whitelist.bed input.bam
```

where `mm9-whitelist.bed` contains not blacklisted genomic regions and `input.bam` is the alignment file obtained with `bwa mem`. Resulting filtered alignment files were sorted by coordinate using `samtools sort` and files corresponding to different sequencing lanes of the same library were merged using `samtools merge`.

Sequencing duplicates were marked using `picard`:

```bash
java -Xmx3G -jar picard.jar MarkDuplicates I=input.bam O=output.bam M=log.markdup.txt
```

were `input.bam` is the merged and sorted alignment file, with duplicates marked as `output.bam` and processing log details written to `log.markdup.txt`. Resulting files were indexed using `samtools index`.

Details about the insert size were obtained using `picard` too:

```bash
java -Xmx2G -jar picard.jar CollectInsertSizeMetrics I=input.bam O=InsertSize.txt H=InsertSize.pdf AS=true VALIDATION_STRINGENCY=SILENT
```


## Nucleosome detection

### [iNPS](http://www.picb.ac.cn/hanlab/iNPS.html)

iNPS requires `.bed` files as input. First convert processed `.bam` files above into `.bed` fortmat. iNPS requires python3.


### [danpos]

Under construction ...

- https://github.com/sblab-bioinformatics/projects/blob/master/20151125_5fC_nucleosome/20161122_mnase/20161122_mnase.md
- https://github.com/sblab-bioinformatics/projects/blob/master/20151125_5fC_nucleosome/20170516_heart_bw_iNPS/20170516_heart_bw_iNPS.md
- http://10.20.192.25:3000/projects/20150309-gordon-redbs-mesc/wiki/14-08-2015_MNase-Seq_heart_and_hindbrain_SLX-10219



## TODO




