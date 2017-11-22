# Analysis of sequencing of genomic DNA fed with only C, 5mC, 5hmC or 5fC.

Each sequencing run is done by duplicate, plus an input library without treatment.

## Trim adapters, align, filter and merge lanes

Trim adapters

```bash
#! /usr/bin/env bash
declare -A new_names=( \
        ["gDNA_C"]="gDNA_C" \
        ["gDNA_mC"]="gDNA_mC" \
        ["gDNA_hmC"]="gDNA_hmC" \
        ["gDNA_fC"]="gDNA_fC" \
        ["gDNA_C_2"]="gDNA_C_2" \
        ["gDNA_mC_2"]="gDNA_mC_2" \
        ["gDNA_hmC_2"]="gDNA_hmC_2" \
        ["gDNA_fC_2"]="gDNA_fC_2" \
        ["gDNA_C_input"]="gDNA_C_input" \
        ["gDNA_mC_input"]="gDNA_mC_input" \
        ["gDNA_hmC_input"]="gDNA_hmC_input" \
        ["gDNA_fC_input"]="gDNA_fC_input" \
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

Align,

```bash
#! /usr/bin/env bash
declare -A new_names=( \
        ["gDNA_C"]="gDNA_C" \
        ["gDNA_mC"]="gDNA_mC" \
        ["gDNA_hmC"]="gDNA_hmC" \
        ["gDNA_fC"]="gDNA_fC" \
        ["gDNA_C_2"]="gDNA_C_2" \
        ["gDNA_mC_2"]="gDNA_mC_2" \
        ["gDNA_hmC_2"]="gDNA_hmC_2" \
        ["gDNA_fC_2"]="gDNA_fC_2" \
        ["gDNA_C_input"]="gDNA_C_input" \
        ["gDNA_mC_input"]="gDNA_mC_input" \
        ["gDNA_hmC_input"]="gDNA_hmC_input" \
        ["gDNA_fC_input"]="gDNA_fC_input" \
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
        | samtools sort -@8 - -o ${bam_folder}/${rr}.bam &&\
        samtools index ${bam_folder}/${rr}.bam"

    echo ${cmd}
    sbatch -o align_${rr}.log -J al_${inp_folder}_${rr} --wrap "${cmd}"
    done

done
```

Merge lanes, filter, remove duplicates.

```bash
#! /usr/bin/env bash
declare -A new_names=( \
        ["gDNA_C"]="gDNA_C" \
        ["gDNA_mC"]="gDNA_mC" \
        ["gDNA_hmC"]="gDNA_hmC" \
        ["gDNA_fC"]="gDNA_fC" \
        ["gDNA_C_2"]="gDNA_C_2" \
        ["gDNA_mC_2"]="gDNA_mC_2" \
        ["gDNA_hmC_2"]="gDNA_hmC_2" \
        ["gDNA_fC_2"]="gDNA_fC_2" \
        ["gDNA_C_input"]="gDNA_C_input" \
        ["gDNA_mC_input"]="gDNA_mC_input" \
        ["gDNA_hmC_input"]="gDNA_hmC_input" \
        ["gDNA_fC_input"]="gDNA_fC_input" \
)

whitelist="mm9-whitelist.bed"
picard_exe="~/programs/picard-2.8.2.jar"

for folder in "${!new_names[@]}"
do
    inp_folder="bam_"${folder}
    out_folder="bam_"${folder}
    who=( ${new_names[$folder]}_L???.bam )
    pairs=()
    for ff in ${inp_folder}/${who[@]}
    do
        pairs+=("${ff}")
    done

    name=${new_names[$folder]}
    cmd="samtools merge -f -@ 20 ${inp_folder}/merged_${name}.tmp.bam ${pairs[@]} && \
        java -Xmx3G -jar ${picard_exe} MarkDuplicates VALIDATION_STRINGENCY=SILENT I=${inp_folder}/merged_${name}.tmp.bam O=${inp_folder}/merged_${name}.bam M=${inp_folder}/${name}.markdup.txt &&\
        echo ' Now let us remove  ' && \
        samtools view ${inp_folder}/merged_${name}.bam -f 3 -F 3840 -q 10 -b -L ${whitelist}  \
        | samtools sort -@ 20 -T /scratcha/sblab/portel01/tmp/${name} -o ${inp_folder}/${name}.bam - && \
        echo ' Done sorting  ' && \
        samtools index ${inp_folder}/${name}.bam && \
        echo ' Done with all ' && \
        rm -fr ${inp_folder}/merged_${name}.bam && echo 'Done' "

    #echo ${cmd}
    sbatch -o merge_${name}.log -J m_${inp_folder}_${name} --wrap "${cmd}"
done
```

### Generate `BigWig` files for later usage (e.g. coverage around 5fC sites)

We use the `extendReads` options from deeptools' `bamCoverage` to the the coverage across the read. We have previously removed duplicates and made sure that all reads are properly paired.

```bash
#! /usr/bin/env bash
declare -A new_names=( \
        ["gDNA_C"]="gDNA_C" \
        ["gDNA_mC"]="gDNA_mC" \
        ["gDNA_hmC"]="gDNA_hmC" \
        ["gDNA_fC"]="gDNA_fC" \
        ["gDNA_C_2"]="gDNA_C_2" \
        ["gDNA_mC_2"]="gDNA_mC_2" \
        ["gDNA_hmC_2"]="gDNA_hmC_2" \
        ["gDNA_fC_2"]="gDNA_fC_2" \
        ["gDNA_C_input"]="gDNA_C_input" \
        ["gDNA_mC_input"]="gDNA_mC_input" \
        ["gDNA_hmC_input"]="gDNA_hmC_input" \
        ["gDNA_fC_input"]="gDNA_fC_input" \
)

for folder in "${!new_names[@]}"
do
    out_folder="bw_"${folder}
    inp_folder="bam_"${folder}

    if [ ! -d "${out_folder}" ]; then
        mkdir ${out_folder}
    fi

    cmd="bamCoverage -b ${inp_folder}/${folder}.bam -o ${out_folder}/${folder}.bw -of bigwig --binSize 1 -p 20 --normalizeUsingRPKM --extendReads --ignoreForNormalization chrX chrM chrY"
    echo ${cmd}
    sbatch -o bw_${folder}.log -J bw_${folder} --wrap "${cmd}"

done
```
