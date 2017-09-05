
## Software

- [fastqc v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [bwa v0.7.12-r1044](http://bio-bwa.sourceforge.net/)
- [samtools v1.2](http://samtools.sourceforge.net/)
- [picard v1.138](https://broadinstitute.github.io/picard/)
- [python v3.5](https://www.python.org/)
- [bedtools v2.25.0](http://bedtools.readthedocs.io/en/latest/)
- common unix tools (e.g. awk, sort and uniq)
- [iNPS_v1.2.2](http://www.picb.ac.cn/hanlab/iNPS.html)


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

iNPS requires `python3` and also that processed `.bam` files above are converted into `.bed` format for input. First, duplicate reads in processed `.bam` files were removed using `samtools view`. Second, sort alignments by read name with `samtools sort` and fill in mate coordinates and related flags using `samtools fixmate`. Third, convert resulting `.bam` files into `.bed` with `bamToBed` and process output for appropiate start and end coordinates depending on read and mate.

```bash
samtools view -u -f 2 -F 1024 input_processed.bam | \
samtools sort -T $outinps -O bam -@8 -n - | \
samtools fixmate -O bam - - | \
bamToBed -bedpe | \
awk -v OFS='\t' '{if(\$2 < \$5){start=\$2} else {start=\$5} if(\$3 > \$6){end=\$3} else {end=\$6} print \$1, start, end}' | \
sort -k1,1 -k2,2n -k3,3n > output.bed"
```

Running iNPS:

```bash
python3.5 iNPS_V1.2.2.py -i input.bed -o ../iNPS/input --s_p=p
```

where `input.bed` is the output file obtained above and `--s_p=p` indicates paired end sequencing data and `--pe_max=230`.


### [danpos](https://sites.google.com/site/danposdoc/)

Under construction ...

