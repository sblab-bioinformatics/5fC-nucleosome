

# ENCODE dataset analyzed containing histone marks 

| Tissue | Antibody | UCSC Accession | 
| ----- | ---- | ---- | 
| Hindbrain | H3K27ac | ENCFF203QTV |
| Hindbrain | H3K4me1 | ENCFF542GAS |
| Heart | H3K27ac | ENCFF954URD | 
| Heart | H3K4me1 | ENCFF737FNO | 
| 

## CpG islands for mm9

We downloaded the CpG islands from [here](https://genome.ucsc.edu/cgi-bin/hgTables) setting "group:All Tracks" and

- "track:CpG Islands" and "table:cpgIslandExt" (mm9.cpgIslandExt.bed)
or
- "track:Unmasked CpG" and "table:cpgIslandExtUnmasked" (mm9.cpgIslandExtUnmasked.bed)

see difference between `cpgIslandExt` and `cpgIslandExtUnmasked` by clicking on `describe table schema` then downloaded the file locally.

```bash
wc -l mm9.cpgIslandExt.bed # 16026
wc -l mm9.cpgIslandExtUnmasked.bed # 23686
```
