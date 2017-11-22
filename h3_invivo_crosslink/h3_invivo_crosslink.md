# Chip-seq H3 histone crosslinked to 5fC

## Trim adapters, align and merge lanes

Trim the adapters,

```bash
#! /usr/bin/env bash

declare -A new_names=( \
        ["h3_reduced_1"]="h3_reduced_1" \
        ["h3_reduced_2"]="h3_reduced_2" \
)

for folder in "${!new_names[@]}"
do
    pairs=()
    who=( ${new_names[$folder]}_L???_??_???.fastq.gz )
    for ff in ${folder}/${who[@]}
    do
        pairs+=("${ff%_??_???.*}")
    done
    sorted_unique_ids=($(echo "${pairs[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
    out_folder="trim_"${folder}

    if [ ! -d "${out_folder}" ]; then
        mkdir ${out_folder}
    fi

    for rr in  "${sorted_unique_ids[@]}"
    do
        cmd="trim_galore -q 10 --stringency 8 -o ${out_folder} --paired ${rr}_R1_001.fastq.gz ${rr}_R2_001.fastq.gz"
        echo ${cmd}
        sbatch -o ${rr}.log -J ${inp_folder}_${rr} --wrap "${cmd}"
    done
done
exit
```

Align all the reads,

```bash
#! /usr/bin/env bash


# let's have some fun with associative arrays
declare -A new_names=( \
        ["h3_reduced_1"]="h3_reduced_1" \
        ["h3_reduced_2"]="h3_reduced_2" \
)

ref="~/genomes/mus_musculus/mm9/Sequence/BWAIndex/genome.fa"
whitelist="mm9-whitelist.bed"
picard_exe="~/programs/picard-2.8.2.jar"

for folder in "${!new_names[@]}"
do
    inp_folder="trim_"${folder}
    bam_folder="bam_"${folder}
    who=( ${new_names[$folder]}_L???_??_???_val_?.fq.gz )
    pairs=()
    for ff in ${inp_folder}/${who[@]}
    do
        filename=$(basename "${ff}")
        short=${filename%_??_???_val_?.fq.gz}
        pairs+=("${short}")
    done
    sorted_unique_ids=($(echo "${pairs[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
    out_folder="bam_"${folder}

    if [ ! -d "${out_folder}" ]; then
        mkdir ${out_folder}
    fi

    for rr in  "${sorted_unique_ids[@]}"
    do
    cmd="bwa mem -M -t 8 $ref ${inp_folder}/${rr}_R1_001_val_1.fq.gz ${inp_folder}/${rr}_R2_001_val_2.fq.gz \
        | samtools sort -@8 - -o ${bam_folder}/${rr}.tmp.bam  &&
        java -Xmx3G -jar ${picard_exe} MarkDuplicates VALIDATION_STRINGENCY=SILENT I=${bam_folder}/${rr}.tmp.bam O=${bam_folder}/${rr}.bam M=${bam_folder}/${rr}.markdup.txt &&
        samtools index ${bam_folder}/${rr}.bam"

    echo ${cmd}
    sbatch -o align_${rr}.log -J al_${inp_folder}_${rr} --wrap "${cmd}"
    done
done
```

Merge lanes, sort and filter by quality and remove duplicates.

```bash
#! /usr/bin/env bash

declare -A new_names=( \
        ["h3_reduced_1"]="h3_reduced_1" \
        ["h3_reduced_2"]="h3_reduced_2" \
)

picard_exe="~/programs/picard-2.8.2.jar"
whitelist="mm9-whitelist.bed"

for folder in "${!new_names[@]}"
do
    inp_folder="bam_"${folder}
    out_folder="bw_"${folder}
    who=( ${new_names[$folder]}_L???.bam )
    pairs=()
    for ff in ${inp_folder}/${who[@]}
    do
        pairs+=("${ff}")
    done

    name=${new_names[$folder]}
    cmd="samtools merge -f -@ 20 ${inp_folder}/merged_${name}.bam ${pairs[@]} && \
        java -Xmx3G -jar ${picard_exe} MarkDuplicates VALIDATION_STRINGENCY=SILENT I=${inp_folder}/merged_${name}.bam O=${inp_folder}/${name}.unf.bam M=${inp_folder}/${name}.markdup.txt  &&\
        samtools view ${inp_folder}/${name}.unf.bam -f 3 -F 3840 -q 10 -b -L ${whitelist}  \
        | samtools sort -@ 20 -T /scratcha/sblab/portel01/tmp/${name} -o ${inp_folder}/${name}.bam - && \
        samtools index ${inp_folder}/${name}.bam && \
        rm -fr ${inp_folder}/merged_${name}.bam ${inp_folder}/${name}.unf.bam"
    echo ${cmd}
    sbatch -o merge_${name}.log -J m_${inp_folder}_${name} --wrap "${cmd}"
done
```

## Peak calling

We use MACS2 for calling peaks, and retain the narrowpeak files as regions with H3 histone crosslinked to 5fC.

```bash
#! /usr/bin/env bash

for rep in 1 2
do

out_folder=peak_calling_macs2_${rep}

if [ ! -d "${out_folder}" ]; then
    mkdir ${out_folder}
fi
cd ${out_folder}

bam_red=bam_h3_reduced_${rep}/h3_reduced_${rep}.bam

cmd="macs2 callpeak --keep-dup all --bdg -f BAMPE -t ../${bam_red}  -g mm -n h3_invivo_no_inp_${rep}  "
sbatch -o mcs2_${rep}.log -J mcs_${rep} --wrap "${cmd}"

cd ../
done
```
