
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
