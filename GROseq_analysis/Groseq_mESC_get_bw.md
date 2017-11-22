# Generate `BigWig` files from Gro-seq reads

The raw data has been taken from [Wang et al., Nature 2015](https://www.ncbi.nlm.nih.gov/pubmed/26123024)

The GEO for the GRO-seq data is GSE64748. There are 5 datasets (with duplicates) for WT and for TDG-KO.

1.  (WT/KO)-DRB3H --> Time zero of GRO-seq
2.  (WT/KO)-WASH10M --> 10 min after removing the blockage of PolII
3.  (WT/KO)-WASH20M --> 20 min after removing the blockage of PolII
4.  (WT/KO)-WASH30M --> 30 min after removing the blockage of PolII
5.  (WT/KO)-NODRB --> Reference for no blockage of PolII

The protocol is partially based on the original paper, updated based on [Nagari et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5522910/)

### Get the data, trimming & alignment

##### Download `.fastq` from GEO

Obtain the `SRA` number from GEO. Then, using `SRAToolkit`, download the data and rename the files.

```bash
#! /usr/bin/env bash
declare -A runs=(
    ["SRR1744498"]="WTNODRB_1"
    ["SRR1744499"]="WTNODRB_2"
    ["SRR1744500"]="WTDRB3H_1"
    ["SRR1744501"]="WTDRB3H_2"
    ["SRR1744502"]="WTWASH10M_1"
    ["SRR1744503"]="WTWASH10M_2"
    ["SRR1744504"]="WTWASH20M_1"
    ["SRR1744505"]="WTWASH20M_2"
    ["SRR1744506"]="WTWASH30M_1"
    ["SRR1744507"]="WTWASH30M_2"
    ["SRR1744508"]="KONODRB_1"
    ["SRR1744509"]="KONODRB_2"
    ["SRR1744510"]="KODRB3H_1"
    ["SRR1744511"]="KODRB3H_2"
    ["SRR1744512"]="KOWASH10M_1"
    ["SRR1744513"]="KOWASH10M_2"
    ["SRR1744514"]="KOWASH20M_1"
    ["SRR1744515"]="KOWASH20M_2"
    ["SRR1744516"]="KOWASH30M_1"
    ["SRR1744517"]="KOWASH30M_2"
)

prefetch="~/programs/sratoolkit.2.8.1-3-centos_linux64/bin/prefetch"
fqdump="~/programs/sratoolkit.2.8.1-3-centos_linux64/bin/fastq-dump"

for run_id in ${!runs[@]}
do
    echo ${run_id}
    sbatch -J ${run_id} -o ${run_id}.log --mem 1024 --wrap "
    ${prefetch} -v ${run_id}
    ${fqdump} --outdir fastq --gzip  --split-files ~/ncbi/public/sra/$run_id.sra &&
    rm -fr ~/ncbi/public/sra/$run_id.sra "
done

# Rename files
for run_id in ${!runs[@]}
do
    echo ${run_id}, ${runs[$run_id]}
    mv fastq/${run_id}_1.fastq.gz fastq/${runs[$run_id]}.fastq.gz
done
```

##### Remove adapters and poly-A tails

```bash
#! /usr/bin/env bash

declare -A new_names=(
    ["WTNODRB_1"]="WTNODRB",
    ["WTNODRB_2"]="WTNODRB",
    ["WTDRB3H_1"]="WTDRB3H",
    ["WTDRB3H_2"]="WTDRB3H"
    ["WTWASH10M_1"]="WTWASH10M"
    ["WTWASH10M_2"]="WTWASH10M"
    ["WTWASH20M_1"]="WTWASH20M"
    ["WTWASH20M_2"]="WTWASH20M"
    ["WTWASH30M_1"]="WTWASH30M"
    ["WTWASH30M_2"]="WTWASH30M"
    ["KONODRB_1"]="KONODRB"
    ["KONODRB_2"]="KONODRB"
    ["KODRB3H_1"]="KODRB3H"
    ["KODRB3H_2"]="KODRB3H"
    ["KOWASH10M_1"]="KOWASH10M"
    ["KOWASH10M_2"]="KOWASH10M"
    ["KOWASH20M_1"]="KOWASH20M"
    ["KOWASH20M_2"]="KOWASH20M"
    ["KOWASH30M_1"]="KOWASH30M"
    ["KOWASH30M_2"]="KOWASH30M"
)

for folder in "${!new_names[@]}"
do
    out_folder="new_trim_"${folder}
    if [ ! -d "${out_folder}" ]; then
        mkdir ${out_folder}
    fi
    cmd="cutadapt -a "TCGTATGCCGTCTTCTGCTTG" -q 20 -e 0.2 -o ${out_folder}/${folder}.fastq.gz  fastq/${folder}.fastq.gz &&
    cutadapt -a "AAAAAAAAAAAAAAAAAAAA" -z -e 0.10 -m 16 -o ${out_folder}/${folder}_nopolyA.fastq.gz  ${out_folder}/${folder}.fastq.gz && rm ${out_folder}/${folder}.fastq.gz"
    sbatch -o ${folder}_1.log -J ${out_folder} --wrap "${cmd}"
done
```

##### Align to reference genome (mm9)

We are using `bwa`, but we obtain very similar results using `bowtie2`.
After alignement, we sort and index the `bam` files.

```bash
#! /usr/bin/env bash
declare -A new_names=(
    ["WTNODRB_1"]="WTNODRB",
    ["WTNODRB_2"]="WTNODRB",
    ["WTDRB3H_1"]="WTDRB3H",
    ["WTDRB3H_2"]="WTDRB3H"
    ["WTWASH10M_1"]="WTWASH10M"
    ["WTWASH10M_2"]="WTWASH10M"
    ["WTWASH20M_1"]="WTWASH20M"
    ["WTWASH20M_2"]="WTWASH20M"
    ["WTWASH30M_1"]="WTWASH30M"
    ["WTWASH30M_2"]="WTWASH30M"
    ["KONODRB_1"]="KONODRB"
    ["KONODRB_2"]="KONODRB"
    ["KODRB3H_1"]="KODRB3H"
    ["KODRB3H_2"]="KODRB3H"
    ["KOWASH10M_1"]="KOWASH10M"
    ["KOWASH10M_2"]="KOWASH10M"
    ["KOWASH20M_1"]="KOWASH20M"
    ["KOWASH20M_2"]="KOWASH20M"
    ["KOWASH30M_1"]="KOWASH30M"
    ["KOWASH30M_2"]="KOWASH30M"
)

ref="~/reference_data/genomes/mus_musculus/mm9/Sequence/BWAIndex/genome.fa"

for folder in "${!new_names[@]}"
do
    out_folder="bams_"${folder}
    inp_folder="trim_"${folder}

    if [ ! -d "${out_folder}" ]; then
        mkdir ${out_folder}
    fi

    cmd=" bwa aln -n 2 -l 25 -t 8 ${ref} ${inp_folder}/${folder}_nopolyA.fastq.gz > ${out_folder}/${folder}.sai &&
    bwa samse ${ref} -n 1 ${out_folder}/${folder}.sai ${inp_folder}/${folder}_nopolyA.fastq.gz > ${out_folder}/${folder}.sam &&
    samtools view -bh -S ${out_folder}/${folder}.sam > ${out_folder}/${folder}.unsorted.bam &&
    samtools sort ${out_folder}/${folder}.unsorted.bam  ${out_folder}/${folder}.bam &&
    samtools index ${out_folder}/${folder}.bam &&
    rm -fr ${out_folder}/${folder}.sai  ${out_folder}/${folder}.sam ${out_folder}/${folder}.unsorted.bam "
    sbatch -o ${folder}.log -J ${out_folder} --wrap "${cmd}"

done
```

##### Merge technical duplicates and remove 'optical' duplicates

We are merging the tehcnical duplicates as it was done in the `Wang et al.` paper on the basis of good correlation between the two sets.

```bash
#! /usr/bin/env bash

samples=( WTNODRB WTDRB3H WTWASH10M WTWASH20M WTWASH30M KONODRB KODRB3H KOWASH10M KOWASH20M KOWASH30M )

picard_exe="~/programs/picard-2.8.2/picard-2.8.2.jar"

for s in ${samples[@]}
do
    out_folder="bams_merged_"${s}

    if [ ! -d "${out_folder}" ]; then
        mkdir ${out_folder}
    fi
    cmd=" samtools merge -f -@ 20 ${out_folder}/${s}.unsorted.bam bams_${s}_1/${s}_1.bam bams_${s}_2/${s}_2.bam &&
          samtools sort -@ 20 -T /scratcha/sblab/portel01/tmp/${s} -o ${out_folder}/${s}.tmp.bam ${out_folder}/${s}.unsorted.bam  && \
    java -Xmx3G -jar ${picard_exe} MarkDuplicates VALIDATION_STRINGENCY=SILENT I=${out_folder}/${s}.tmp.bam O=${out_folder}/${s}.bam M=${out_folder}/${s}.markdup.txt REMOVE_DUPLICATES=true &&\
        samtools index ${out_folder}/${s}.bam && \
        rm -fr ${out_folder}/${s}.unsorted.bam "
    echo ${cmd}
    sbatch -o merge_${s}.log -J m_${s} --wrap "${cmd}"
done
```

##### Obtain `BigWig` files for coverage Analysis

We create on file for the forward and one file for the reverse strand coveraged. In this particular case, `Deeptools'` `bamCoverage` tool
does the opposite of what would be expected, so when you tell it to filter the forward strand, this actually produces the reverse one.

```bash
#! /usr/bin/env bash

samples=( WTNODRB WTDRB3H WTWASH10M WTWASH20M WTWASH30M KONODRB KODRB3H KOWASH10M KOWASH20M KOWASH30M )

for s in ${samples[@]}
do
    out_folder="bw_"${s}
    if [ ! -d "${out_folder}" ]; then
        mkdir ${out_folder}
    fi
    ## BE CAREFULL HERE, BAMCOVERAGE MAPS THE STRAND INFORMATION IN THE OTHER STRAND FOR THIS
    ## PARTICULAR CASE
    sbatch -J ${s} -o "bam2bw_"${s} --wrap "bamCoverage -b bams_merged_${s}/${s}.bam -o ${out_folder}/${s}_r.bw -of bigwig --binSize 1 -p 8 --normalizeUsingRPKM --filterRNAstrand forward"
    sbatch -J ${s} -o "bam2bw_"${s} --wrap "bamCoverage -b bams_merged_${s}/${s}.bam -o ${out_folder}/${s}_f.bw -of bigwig --binSize 1 -p 8 --normalizeUsingRPKM --filterRNAstrand reverse"
done
```
