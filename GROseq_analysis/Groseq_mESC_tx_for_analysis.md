# Select active transcripts for metagene analysis

From the aligned `bams`, we calculated the RPKMs on mm9 transcripts using this [R script](groseq_rpkm.R)

```R
source("http://bioconductor.org/biocLite.R")
biocLite("org.Mm.eg.db")
biocLite("TxDb.Mmusculus.UCSC.mm9.knownGene")
library(groHMM)
library(GenomicFeatures)
library(edgeR)
library(org.Mm.eg.db)
options(mc.cores=getCores(10))
# WT
Reads_wtnodrb_1 <- as(readGAlignments("../bams_WTNODRB_1/WTNODRB_1.bam"), "GRanges")

library(TxDb.Mmusculus.UCSC.mm9.knownGene)
kgdb <- TxDb.Mmusculus.UCSC.mm9.knownGene


# For known genes,
kg <- transcripts(kgdb, columns=c("gene_id", "tx_id", "tx_name"))

txlen <- na.omit(transcriptLengths(kgdb))

map <- select(org.Mm.eg.db,
              keys=unique(unlist(mcols(kg)$gene_id)),
              columns=c("SYMBOL"), keytype=c("ENTREZID"))
missing <- elementNROWS(mcols(kg)[,"gene_id"]) == 0
kg <- kg[!missing,]
inx <- match(unlist(mcols(kg)$gene_id), map$ENTREZID)
mcols(kg)$symbol <- map[inx,"SYMBOL"]
kgLimit <- limitToXkb(kg)
ctWTNODRB_1  <- countOverlaps(kgLimit, Reads_wtnodrb_1)


dfnames <- data.frame( cbind(kg$tx_name, kg$symbol, start(ranges(kg)), end(ranges(kg)) ))
dfnames <- data.frame( cbind(kg$tx_name, kg$symbol, as.character(seqnames(kg)), start(ranges(kg)), end(ranges(kg)),  as.character(strand(kg)) ))
colnames(dfnames) <- c("tx_name", "gene_name", "chr", "start", "end", "strand")

rpkm.d <- rpkm(ctWTNODRB_1, txlen$tx_len)
rpkm.d<- cbind(kg$tx_name, kg$symbol, as.character(seqnames(kg)), start(ranges(kg)), end(ranges(kg)),
                              as.character(strand(kg)), rpkm.d)
colnames(rpkm.d) <- c("tx_name", "gene_name", "chr", "start", "end", "strand", "rpkm")
write.table(rpkm.d, 'rpkm_WTNODRB.txt', row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
```

which could be run e.g.

```bash
Rscript groseq_rpkm.R
```

As in the original Wang et al. paper, for the `metagene` analysis we only keep those transcripts which show more the 0.5 RPKMs in the `WTNOBRB` sample.

```bash
egrep -v "tx_name"  rpkm_WTNODRB.txt | awk -v "OFS=\t" '($NF > 0.5) {print $3, $4, $5, $2, "0", $6}' | sortBed -i -  > tx_wtnodrb_lg_05rpkm.bed
```

Prepare `bed` files denoting regions of interest for our metagene analysis. The 5fC and 5caC sites were obtained from [Song et al. Cell 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3657391/)

The `h3_invivo_1_peaks_*.narroPeak` files were obtained from our Chip-seq of H3 crosslinked to 5fC, the protocol is described [here](../h3_invivo_crosslink/h3_invivo_crosslink.md).

```bash
#! /usr/bin/env bash

# with rpkm>0.5
bed=tx_wtnodrb_lg_05rpkm.bed
name_bed="${bed%.*}"
fC_mESC=5fC_mESC_kd.bed
caC_mESC=5caC_mESC_kd.bed

# union H3-Chip-seq 1_2
cat h3_invivo_1_peaks.narrowPeak h3_invivo_2_peaks.narrowPeak  > union_h3_1_2.bed
bedtools intersect -a ${bed} -b union_h3_1_2.bed -u | sortBed -i - > ${name_bed}_union_h3_1_3.bed
# union H3-Chip-seq 1_2  that also overlap with known 5fC sites
bedtools intersect -a union_h3_1_2.bed -b ${fC_mESC}  -u | sortBed -i - > union_h3_1_2_overlap5fC.bed
bedtools intersect -a ${bed} -b union_h3_1_2_overlap5fC.bed -u | sortBed -i - > ${name_bed}_union_h3_1_2_overlap5fC.bed

# Sites unique to 5fC, or unique to 5caC or without 5fC
bedtools intersect -a ${bed} -b ${fC_mESC}  -u   | sortBed -i -  > ${name_bed}_5fC.bed
bedtools intersect -a ${bed} -b ${caC_mESC} -u   | sortBed -i -  > ${name_bed}_5caC.bed
bedtools intersect -a ${bed} -b ${fC_mESC}  -v -wa   | sortBed -i -  > ${name_bed}_no5fC.bed
# concat fC and caC to get sites without 5caC and without 5fC
cat  ${fC_mESC} ${caC_mESC} > 5fC_and_5caC_mES_kd.bed
bedtools intersect -a ${bed} -b 5fC_and_5caC_mES_kd.bed -v -wa   | sortBed -i -  > ${name_bed}_no5fC_no5caC.bed

# 5fC without 5caC, 5caC without 5fC
bedtools intersect -a ${fC_mESC} -b ${caC_mESC} -v -wa > 5fC_without_5caC_mESC.bed
bedtools intersect -a ${caC_mESC} -b ${fC_mESC} -v -wa > 5caC_without_5fC_mESC.bed

bedtools intersect -a ${bed} -b 5fC_without_5caC_mESC.bed -u   | sortBed -i -  > ${name_bed}_5fC_no5caC.bed
bedtools intersect -a ${bed} -b 5caC_without_5fC_mESC.bed -v -wa   | sortBed -i -  > ${name_bed}_5caC_no5fC.bed

# Now filter per strand
declare -a proces=( ${name_bed}_5fC.bed ${name_bed}_5caC.bed ${name_bed}_no5fC.bed ${name_bed}_no5fC_no5caC.bed  ${bed}
    ${name_bed}_5fC_no5caC.bed ${name_bed}_5caC_no5fC.bed ${name_bed}_union_h3_1_2_overlap5fC.bed  ${name_bed}_union_h3_1_2.bed )
for f in ${proces[@]}
do
    name=${f%.*}
    cat ${name}.bed | awk '($6 == "+")' > ${name}_f.bed
    cat ${name}.bed | awk '($6 == "-")' > ${name}_r.bed
done
exit
```
