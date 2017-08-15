# Obtain 5fC sites in heart and hindbrain, both for WT and TDG-KO

Prepare a set of `bed` files containing the genomic regions for 5fC in
heart and hindbrain, both for WT and for TDG-KO.

Create sets of 5fC sites that are a consensus between replicates, as well as
common or unique between tissues/conditions.

### Hindbrain

```bash
#! /usr/bin/env bash

# WT
# sort, merge and count the WT
zcat  \
grm069_TDG_WT_O3_hindbrain.clean_peaks.narrowPeak.gz \
grm070_TDG_WT_Q3_hindbrain.clean_peaks.narrowPeak.gz \
| sort -k1,1 -k2,2n -s \
| mergeBed  -c 1 -o count | gzip  > 5fC_merged_WT.bed.gz

# select consenus WT, column $4 == 2 if present in both files
zcat 5fC_merged_WT.bed.gz | awk '($4==2) {print $0}' > 5fC_consensus_WT.bed
gzip 5fC_consensus_WT.bed

# Repeat for KO
zcat \
grm064_TDG_KO_2A2_hindbrain.clean_peaks.narrowPeak.gz \
grm065_TDG_KO_O4_hindbrain.clean_peaks.narrowPeak.gz \
grm066_TDG_KO_Q2_hindbrain.clean_peaks.narrowPeak.gz \
grm067_TDG_KO_Z1_hindbrain.clean_peaks.narrowPeak.gz \
| sort -k1,1 -k2,2n -s | mergeBed  -c 1 -o count | gzip > 5fC_merged_KO.bed.gz

# select consenus KO, column $4 == 4 if present in both files
zcat 5fC_merged_KO.bed.gz | awk '($4==4) {print $0}' > 5fC_consensus_KO.bed
# wc -l 5fC_consensus_KO.bed --> 116506
gzip 5fC_consensus_KO.bed

# Now intersect both KO and WT to find regions common in both.
bedtools intersect -a 5fC_consensus_WT.bed.gz -b 5fC_consensus_KO.bed.gz -r -f 0.4 \
| gzip >  5fC_common_WT_KO.bed.gz

# unique for WT
bedtools intersect  -a 5fC_consensus_WT.bed.gz -b \
grm065_TDG_KO_O4_hindbrain.clean_peaks.narrowPeak.gz \
grm064_TDG_KO_2A2_hindbrain.clean_peaks.narrowPeak.gz \
grm066_TDG_KO_Q2_hindbrain.clean_peaks.narrowPeak.gz \
grm067_TDG_KO_Z1_hindbrain.clean_peaks.narrowPeak.gz \
-v | gzip > 5fC_sites_unique_to_WT.bed.gz

# unique for KO
bedtools intersect -a 5fC_consensus_KO.bed.gz -b \
grm069_TDG_WT_O3_hindbrain.clean_peaks.narrowPeak.gz \
grm070_TDG_WT_Q3_hindbrain.clean_peaks.narrowPeak.gz \
-v | gzip > 5fC_sites_unique_to_KO.bed.gz

# same number of unique WT but without 5fC
bedtools shuffle -i 5fC_sites_unique_to_WT.bed.gz -excl 5fC_merged_WT.bed.gz \
-g mm9.genome.fa.fai.genomefile | bedtools sort -i - | gzip > non5fC_WT.bed.gz
# same number of unique KO but without 5fC
bedtools shuffle -i 5fC_sites_unique_to_KO.bed.gz -excl 5fC_merged_KO.bed.gz \
-g mm9.genome.fa.fai.genomefile | bedtools sort -i - | gzip > non5fC_KO.bed.gz

# and two more
bedtools shuffle -i 5fC_consensus_WT.bed.gz -excl 5fC_sites_merged_wt_and_ko.bed.gz -g mm9.genome.fa.fai.genomefile | bedtools sort -i - | gzip > non5fC_anywhere_same_size_consensus_WT.bed.gz

bedtools shuffle -i 5fC_consensus_KO.bed.gz -excl 5fC_sites_merged_wt_and_ko.bed.gz -g mm9.genome.fa.fai.genomefile | bedtools sort -i - | gzip > non5fC_anywhere_same_size_consensus_KO.bed.gz
```

### Heart

```bash
#! /usr/bin/env bash

WT_1="GSM2052020_heart_WT_i9.bed.gz"
WT_2="GSM2052018_heart_WT_i6.bed.gz"
KO_1="GSM2052016_heart_KO_i5.bed.gz"
KO_2="GSM2052014_heart_KO_i1.bed.gz"

# WT
# sort, merge and count the WT
zcat  \
${WT_1} \
${WT_2} \
| sort -k1,1 -k2,2n -s \
| mergeBed  -c 1 -o count | gzip  > 5fC_heart_merged_WT.bed.gz

# select consenus WT, column $4 == 2 if present in both files
zcat 5fC_heart_merged_WT.bed.gz | awk '($4==2) {print $0}' > 5fC_heart_consensus_WT.bed
# wc -l 5fC_heart_consensus_WT.bed  -->  7118
gzip 5fC_heart_consensus_WT.bed

# Repeat for KO

zcat \
${KO_1} \
${KO_2} \
| sort -k1,1 -k2,2n -s | mergeBed  -c 1 -o count | gzip > 5fC_heart_merged_KO.bed.gz

# select consenus KO, column $4 == 2 if present in both files
zcat 5fC_heart_merged_KO.bed.gz | awk '($4==2) {print $0}' > 5fC_heart_consensus_KO.bed
gzip 5fC_heart_consensus_KO.bed

#Now intersect both `KO` and `WT` to find regions common in both.
bedtools intersect -a 5fC_heart_consensus_WT.bed.gz -b 5fC_heart_consensus_KO.bed.gz -r -f 0.4 \
| gzip >  5fC_heart_common_WT_KO.bed.gz

# unique for WT
bedtools intersect  -a 5fC_heart_consensus_WT.bed.gz -b \
${KO_1} \
${KO_2} \
-v | gzip > 5fC_heart_sites_unique_to_WT.bed.gz

# unique for KO
bedtools intersect -a 5fC_heart_consensus_KO.bed.gz -b \
${WT_1} \
${WT_2} \
-v | gzip > 5fC_heart_sites_unique_to_KO.bed.gz

# same number of unique WT but without 5fC_heart
bedtools shuffle -i 5fC_heart_sites_unique_to_WT.bed.gz -excl 5fC_heart_merged_WT.bed.gz \
-g mm9.genome.fa.fai.genomefile | bedtools sort -i - | gzip > non5fC_heart_WT.bed.gz
# same number of unique KO but without 5fC_heart
bedtools shuffle -i 5fC_heart_sites_unique_to_KO.bed.gz -excl 5fC_heart_merged_KO.bed.gz \
-g mm9.genome.fa.fai.genomefile | bedtools sort -i - | gzip > non5fC_heart_KO.bed.gz

zcat ${WT_1} ${WT_2} ${KO_1} ${KO_2} | bedtools sort -i - | bedtools merge -i - | gzip > 5fC_heart_sites_merged_wt_and_ko.bed.gz

# and two more
bedtools shuffle -i 5fC_heart_consensus_WT.bed.gz -excl 5fC_heart_sites_merged_wt_and_ko.bed.gz -g mm9.genome.fa.fai.genomefile | bedtools sort -i - | gzip > non5fC_heart_anywhere_same_size_consensus_WT.bed.gz

bedtools shuffle -i 5fC_heart_consensus_KO.bed.gz -excl 5fC_heart_sites_merged_wt_and_ko.bed.gz -g mm9.genome.fa.fai.genomefile | bedtools sort -i - | gzip > non5fC_heart_anywhere_same_size_consensus_KO.bed.gz
```
