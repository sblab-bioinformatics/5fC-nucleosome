
# Calculating nucleosome coverage (as in RPKM) around a given genomic feature

First, we intersect the location of histone marks / CpG islands with regions 
containing 5fC, as well as those without.  

### Intersect histone marks / CpG islands 

Notice that the paths are particular to the machine.

```bash
#! /usr/bin/env bash

raw_data="/scratcha/sblab/portel01/project/euni_5fc/chipseq_5fC_data/"

wt_5fC=${raw_data}/5fC_consensus_WT.bed.gz
ko_5fC=${raw_data}/5fC_consensus_KO.bed.gz
wt_ko_5fC=${raw_data}/5fC_common_WT_KO.bed.gz

# add chr in front, required becase the annotations have it
chr_wt_5fC=chr_5fC_consensus_WT.bed.gz
chr_ko_5fC=chr_5fC_consensus_KO.bed.gz

chr_wt_ko_5fC=chr_5fC_common_WT_KO.bed.gz

zcat ${wt_5fC} | awk '{print "chr"$0}' | gzip > ${chr_wt_5fC}
zcat ${ko_5fC} | awk '{print "chr"$0}' | gzip > ${chr_ko_5fC}
zcat ${wt_ko_5fC} | awk '{print "chr"$0}' | gzip > ${chr_wt_ko_5fC}

# These are the ref names for the following histone modifications
# already after the lift-over to mm9
H3K27ac=ENCFF203QTV.mm9.bed
H3K4me1=ENCFF542GAS.mm9.bed
# from https://genome.ucsc.edu/cgi-bin/hgTables,
# see "Histone mark and CpG island datasets used" entry in this repository 
CpGi=mm9.cpgIslandExtUnmasked.bed

for modif in H3K27ac H3K4me1 CpGi
do

    # not very elegant, but works
    for method in overlap reverse
    do

    # using the variable as another variable is frowned upon, too bad
    working_bed=${!modif}

    if [ ${method} == reverse ]
    then

        # WT
        bedtools intersect  -v -b ${chr_wt_5fC}  -a ${working_bed} \
        | awk -v kk=${modif} '{printf"%s\t%d\t%d\t%s\n",$1,$2,$3,kk}' | gzip \
        > ${modif}_no5fC_consensus_WT.bed.gz

        # KO
        bedtools intersect -v -b ${chr_ko_5fC}  -a ${working_bed}  \
        | awk -v kk=${modif} '{printf"%s\t%d\t%d\t%s\n",$1,$2,$3,kk}' | gzip \
        > ${modif}_no5fC_consensus_KO.bed.gz

        # WT intersect KO
        bedtools intersect  -v -b ${chr_wt_ko_5fC}  -a ${working_bed} \
        | awk -v kk=${modif} '{printf"%s\t%d\t%d\t%s\n",$1,$2,$3,kk}' | gzip \
        > ${modif}_no5fC_common_WT_KO.bed.gz

    else

        # WT
        bedtools intersect -b ${chr_wt_5fC}  -a ${working_bed} -u -wa -f 0.4 -r   \
        | awk -v kk=${modif} '{printf"%s\t%d\t%d\t%s\n",$1,$2,$3,kk}' | gzip \
        > ${modif}_5fC_consensus_WT.bed.gz

        # KO
        bedtools intersect -b ${chr_ko_5fC}  -a ${working_bed} -u -wa -f 0.4 -r  \
        | awk -v kk=${modif} '{printf"%s\t%d\t%d\t%s\n",$1,$2,$3,kk}' | gzip \
        > ${modif}_5fC_consensus_KO.bed.gz

        # WT intersect KO
        bedtools intersect  -b ${chr_wt_ko_5fC}  -a ${working_bed} -u -wa -f 0.4 -r \
        | awk -v kk=${modif} '{printf"%s\t%d\t%d\t%s\n",$1,$2,$3,kk}' | gzip \
        > ${modif}_5fC_common_WT_KO.bed.gz

    fi


    done
done

```

###  MNase-seq read around CpG islands

Using Deeptools with the previously generated `bed` files, and the `bigwig` files 
from [MNase-seq](../MNase-seq/) alignment. 


```bash
computeMatrix reference-point   -R CpGi_5fC_consensus_WT.bed.gz CpGi_no5fC_consensus_WT.bed.gz \
-S ear030_HindbrainA1.150618.ncbi37.bw ear031_HindbrainC1.150618.ncbi37.bw ear032_HindbrainF1.150618.ncbi37.bw \
-out CpGi_WT_referencepoint.mat.gz  \
--referencePoint "center"        \
-b 1000  -a 1000  -bs 10   --skipZeros  --sortRegions "no"   -p "max"  && \

plotProfile --matrixFile CpGi_WT_referencepoint.mat.gz  -out CpGi_WT_referencepoint.pdf  \
--outFileNameData CpGi_WT_referencepoint.dat        \
--plotType "se"  \
--colors red red red         
--samplesLabel WT_A1 WT_C1 WT_F1  \    
--regionsLabel "CpGi" "CpGi_no5fC" \
-y "RPKM"  --perGroup
```

Finally plot `CpGi_WT_referencepoint.dat`, e.g. like [this](scripts/plot_5fc_no5fc_in_one.ipynb)
