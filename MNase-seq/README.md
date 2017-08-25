
## Software

- [fastqc v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [bwa 0.7.12-r1044](http://bio-bwa.sourceforge.net/)
- [samtools v1.2](http://samtools.sourceforge.net/)
- [picard v1.138](https://broadinstitute.github.io/picard/)
- [bedtools v2.25.0](http://bedtools.readthedocs.io/en/latest/)
- [macs2 v2.1.1.20160226](https://github.com/taoliu/MACS)


## Processing

Sequencing quality was explored using `fastqc`. Alignments to the mouse reference genome `mm9` for paired end sequencing data were generated using `bwa mem`:

```bash
bwa mem -M -t 8 $ref $fq1 $fq2
```

where `$ref` is the fasta reference genome file, and `$fq1` and `$fq2` correspond to R1 and R2 sequencing files respectively.

Unmapped, not primary aligned and supplementary reads were filtered from the resulting alignment `.bam` files using `samtools view`. Reads aligning with low mapping quality or aligning to [blacklisted genomic regions](https://sites.google.com/site/anshulkundaje/projects/blacklists) of the `mm9` reference genome were also filtered out.

```bash
samtools view -@8 -S -u -F2820 -q 5 -L mm9-whitelist.bed $bam
```

where `mm9-whitelist.bed` contains not blacklisted genomic regions and `$bam` is the alignment file obtained with `bwa mem`.

