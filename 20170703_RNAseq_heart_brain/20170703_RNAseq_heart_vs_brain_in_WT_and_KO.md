# Compute the DE between heart and hindbrain tissues, in WT and KO.

Join together the counts from all tissues and conditions,

```bash
#! /usr/bin/env bash

cd /scratcha/sblab/portel01/project/euni_5fc/RNA_seq/Heart_vs_brain_htseqcounts
htseq_matrix.py *.htseq.txt -s '\.htseq$' -s '\..*' | grep -v '^__' > htseq_brain_heart_all.matrix
mkdir ../Heart_vs_brain_DE/
cp htseq_brain_heart_all.matrix ../Heart_vs_brain_DE/
```

Run edgeR on the count data. It requires
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
cd /scratcha/sblab/portel01/project/euni_5fc/RNA_seq/Heart_vs_brain_DE

Rscript run_edgeR_b_h_wt.R
Rscript run_edgeR_b_h_ko.R
```

The relevant files are `expr_diff.brain_vs_heart_ko.edgeR.txt` and `expr_diff.brain_vs_heart_wt.edgeR.txt`.
