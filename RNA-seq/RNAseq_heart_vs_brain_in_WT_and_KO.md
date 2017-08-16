# Compute the DE between heart and hindbrain tissues, in WT and KO.

Join together the counts from all tissues and conditions,

```bash
#! /usr/bin/env bash
htseq_matrix.py *.htseq.txt -s '\.htseq$' -s '\..*' | grep -v '^__' > htseq_brain_heart_all.matrix
```

Run edgeR on the count data to obtain the differential expression and RPKMs. It requires
the scripts [intensityFilter.R](scripts/intensityFilter.R), [TPM.R](scripts/TPM.R),
and [FPKM.R](script/FPKM.R).
The `design.txt` including both tissues and TDG conditions.

```bash
library_id tissue tdg
Brain_KO1_TGACCACT brain ko
Brain_KO2_ACAGTGGT brain ko
Brain_KO3_GCCAATGT brain ko
Brain_WT1_ATCACGTT brain wt
Brain_WT2_CGATGTTT brain wt
Brain_WT3_TTAGGCAT brain wt
Heart_KO1_TAGCTTGT heart ko
Heart_KO2_GGCTACAG heart ko
Heart_KO3_CTTGTACT heart ko
Heart_WT1_CAGATCTG heart wt
Heart_WT2_ACTTGATG heart wt
Heart_WT3_GATCAGCG heart wt
```

```bash
Rscript run_edgeR_b_h_wt.R
Rscript run_edgeR_b_h_ko.R
```

The relevant files are `expr_diff.brain_vs_heart_ko.edgeR.txt` and `expr_diff.brain_vs_heart_wt.edgeR.txt`.
