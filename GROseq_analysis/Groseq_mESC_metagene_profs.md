# Generation of metagene profiles

We will be using `Deeptools` `computeMatrix reference-point` program to generate the raw coverage aligned to the TSS.

For the `fwd` strand

```bash
#! /usr/bin/env bash
if [ "$#" -ne 2 ]; then
    echo "Need a befilename and a folder name"
    exit
fi

bedname=$1
who=$2
out_folder=computeMatrix_${who}

bed=`pwd`/${bedname}

if [ ! -d ${out_folder} ]; then
    mkdir ${out_folder}
fi

# Deeptools
computeMatrix_exe='~/programs/computeMatrix'

std_f="Bigwig_files/bw_WTNODRB/WTNODRB_f.bw"
t0_f="Bigwig_files/bw_WTDRB3H/WTDRB3H_f.bw"
t10_f="Bigwig_files/bw_WTWASH10M/WTWASH10M_f.bw"
t20_f="Bigwig_files/bw_WTWASH20M/WTWASH20M_f.bw"
t30_f="Bigwig_files/bw_WTWASH30M/WTWASH30M_f.bw"

sbatch -J "wt_${who}" -o ${out_folder}/wt_${who}.log --mem 65536 --wrap "${computeMatrix_exe} reference-point \
-R $bed \
-S ${t0_f} ${t10_f} ${t20_f} ${t30_f} ${std_f} \
-out ${out_folder}/metagene_f_wt_${who}.mat.gz \
--referencePoint "TSS" \
-b 20000 \
-a 140000 \
-bs 10 \
--averageTypeBins "mean" \
--nanAfterEnd \
-p "max"
"

ko_std_f="Bigwig_files/bw_KONODRB/KONODRB_f.bw"
ko_t0_f="Bigwig_files/bw_KODRB3H/KODRB3H_f.bw"
ko_t10_f="Bigwig_files/bw_KOWASH10M/KOWASH10M_f.bw"
ko_t20_f="Bigwig_files/bw_KOWASH20M/KOWASH20M_f.bw"
ko_t30_f="Bigwig_files/bw_KOWASH30M/KOWASH30M_f.bw"

##computing matrix
sbatch -J "ko_${who}" -o ${out_folder}/${out_folder}.log --mem 65536 --wrap "${computeMatrix_exe} reference-point \
-R $bed \
-S ${ko_t0_f} ${ko_t10_f} ${ko_t20_f} ${ko_t30_f} ${ko_std_f}  \
-out ${out_folder}/metagene_f_ko_${who}.mat.gz \
--referencePoint "TSS" \
-b 20000 \
-a 140000 \
-bs 10 \
--averageTypeBins "mean" \
--nanAfterEnd \
-p "max"
"
```

and for the `rev` strand, basically the same:

```bash
#! /usr/bin/env bash

if [ "$#" -ne 2 ]; then
    echo "Need a befilename and a folder name"
    exit
fi

bedname=$1
who=$2
out_folder=computeMatrix_${who}

bed=`pwd`/${bedname}

if [ ! -d ${out_folder} ]; then
    mkdir ${out_folder}
fi

# Deeptools
computeMatrix_exe='~/programs/computeMatrix'

std_r="Bigwig_files/bw_WTNODRB/WTNODRB_r.bw"
t0_r="Bigwig_files/bw_WTDRB3H/WTDRB3H_r.bw"
t10_r="Bigwig_files/bw_WTWASH10M/WTWASH10M_r.bw"
t20_r="Bigwig_files/bw_WTWASH20M/WTWASH20M_r.bw"
t30_r="Bigwig_files/bw_WTWASH30M/WTWASH30M_r.bw"


sbatch -J "wt_${who}" -o ${out_folder}/wt_${who}.log --mem 65536 --wrap "${computeMatrix_exe} reference-point \
-R $bed \
-S ${t0_r} ${t10_r} ${t20_r} ${t30_r} ${std_r}  \
-out ${out_folder}/metagene_r_wt_${who}.mat.gz \
--referencePoint "TSS" \
-b 20000 \
-a 140000 \
-bs 10 \
--averageTypeBins "mean" \
--nanAfterEnd \
-p "max"
"

ko_std_r="Bigwig_files/bw_KONODRB/KONODRB_r.bw"
ko_t0_r="Bigwig_files/bw_KODRB3H/KODRB3H_r.bw"
ko_t10_r="Bigwig_files/bw_KOWASH10M/KOWASH10M_r.bw"
ko_t20_r="Bigwig_files/bw_KOWASH20M/KOWASH20M_r.bw"
ko_t30_r="my_Bigwig_files/bw_KOWASH30M/KOWASH30M_r.bw"

##computing matrix
sbatch -J "ko_${who}" -o ${out_folder}/${out_folder}.log --mem 65536 --wrap "${computeMatrix_exe} reference-point \
-R $bed \
-S ${ko_t0_r} ${ko_t10_r} ${ko_t20_r} ${ko_t30_r} ${ko_std_r}  \
-out ${out_folder}/metagene_r_ko_${who}.mat.gz \
--referencePoint "TSS" \
-b 20000 \
-a 140000 \
-bs 10 \
--averageTypeBins "mean" \
--nanAfterEnd \
-p "max"
"
```

Launch all the analysis,

```bash
#! /usr/bin/env bash
bed_name=tx_wtnodrb_lg_05rpkm

for str in f r
do
    ./do_cm_general_${str}.sh ${bed_name}_5caC_no5fC_${str}.bed 5caC_no5fC
    ./do_cm_general_${str}.sh ${bed_name}_no5fC_${str}.bed  no5fC
    ./do_cm_general_${str}.sh ${bed_name}_5fC_no5caC_${str}.bed 5fC_no5caC
    ./do_cm_general_${str}.sh ${bed_name}_5fC_${str}.bed 5fC
    ./do_cm_general_${str}.sh ${bed_name}_${str}.bed all
    ./do_cm_general_${str}.sh ${bed_name}_no5fC_no5caC_${str}.bed no5fC_no5caC
    ./do_cm_general_${str}.sh ${bed_name}_5caC_${str}.bed 5caC
    ./do_cm_general_${str}.sh ${bed_name}_union_h3_1_2_overlap5fC_${str}.bed union_h3_1_2_overlap5fC
    ./do_cm_general_${str}.sh ${bed_name}_union_h3_1_3_${str}.bed union_h3_1_3
echo ${str}
done
```

Once done, we use `plotProfile` to get the averaged metagene profiles. Since we will plot the data ourseflves, we don't chare about the `pdf` it generates, we just need the raw data (--outFileNameData switch).
For the `fwd` strand,

```bash
#! /usr/bin/env bash

if [ "$#" -ne 1 ]; then
    echo "Need a suffix folder name "
    exit
fi

who=$1
out_folder=computeMatrix_${who}


if [ ! -d ${out_folder} ]; then
    echo "Error, computeMatrix folder not present"
    exit
fi

# Deeptools
plotProfile_exe='~/programs/plotProfile'

cd ${out_folder}
sbatch -J "p_wt_${who}" -o plot_wt_${who}.log  --wrap "${plotProfile_exe} \
--matrixFile metagene_f_wt_${who}.mat.gz \
-out metagene_tss_f_${who}.pdf \
--outFileNameData metagene_tss_f_${who}.dat \
--plotType "se" \
--startLabel "start" \
--endLabel "end" \
--samplesLabel t0_f t10_f t20_f t30_f std_f  \
--regionsLabel "Metagene_WT" \
--refPointLabel "Center" \
--legendLocation "upper-right" \
-y "RPKM" \
--perGroup
"

sbatch -J "p_ko_${who}" -o plot_ko_${who}.log --wrap "${plotProfile_exe} \
--matrixFile metagene_f_ko_${who}.mat.gz \
-out metagene_tss_f_ko_${who}.pdf \
--outFileNameData metagene_tss_f_ko_${who}.dat \
--plotType "se" \
--startLabel "start" \
--endLabel "end" \
--samplesLabel t0_f t10_f t20_f t30_f std_f  \
--regionsLabel "Metagene_KO" \
--refPointLabel "Center" \
--legendLocation "upper-right" \
-y "RPKM" \
--perGroup
"
```

and for the `rev`,

```bash
#! /usr/bin/env bash

if [ "$#" -ne 1 ]; then
    echo "Need a suffix folder name "
    exit
fi

who=$1
out_folder=computeMatrix_${who}


if [ ! -d ${out_folder} ]; then
	echo "Error, computeMatrix folder not present"
	exit
fi

# Deeptool
plotProfile_exe='~/programs/plotProfile'

cd ${out_folder}
sbatch -J "p_wt_${who}" -o plot_wt_${who}.log  --wrap "${plotProfile_exe} \
--matrixFile metagene_r_wt_${who}.mat.gz \
-out metagene_tss_r_${who}.pdf \
--outFileNameData metagene_tss_r_${who}.dat \
--plotType "se" \
--startLabel "start" \
--endLabel "end" \
--samplesLabel t0_r t10_r t20_r t30_r std_r \
--regionsLabel "Metagene_WT" \
--refPointLabel "Center" \
--legendLocation "upper-right" \
-y "RPKM" \
--perGroup
"

sbatch -J "p_ko_${who}" -o plot_ko_${who}.log --wrap "${plotProfile_exe} \
--matrixFile metagene_r_ko_${who}.mat.gz \
-out metagene_tss_r_ko_${who}.pdf \
--outFileNameData metagene_tss_r_ko_${who}.dat \
--plotType "se" \
--startLabel "start" \
--endLabel "end" \
--samplesLabel t0_r t10_r t20_r t30_r std_r \
--regionsLabel "Metagene_KO" \
--refPointLabel "Center" \
--legendLocation "upper-right" \
-y "RPKM" \
--perGroup
"
```

Launch the two scripts above,

```bash
#! /usr/bin/env bash
for str in f r
do
    ./do_plot_g_${str}.sh 5caC_no5fC
    ./do_plot_g_${str}.sh no5fC
    ./do_plot_g_${str}.sh 5fC_no5caC
    ./do_plot_g_${str}.sh 5fC
    ./do_plot_g_${str}.sh all
    ./do_plot_g_${str}.sh no5fC_no5caC
    ./do_plot_g_${str}.sh 5caC
    ./do_plot_g_${str}.sh union_h3_1_2_overlap5fC
    ./do_plot_g_${str}.sh union_h3_1_2
echo ${str}
done
```
