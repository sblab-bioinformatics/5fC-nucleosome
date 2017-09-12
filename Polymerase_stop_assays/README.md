
## Software

- [bcl2fastq v2.17.1.14](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html)
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [bowtie2 v2.2.6](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [bwa v0.7.12-r1044](http://bio-bwa.sourceforge.net/)
- [samtools v1.2](http://samtools.sourceforge.net/)
- [python v2.7.10](https://www.python.org/)
- common unix tools (e.g. awk, sort and uniq)
- [R v3.2.1](https://www.r-project.org/)
  - [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
  - [ggplot2](http://ggplot2.org/)
  - [data.table](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html)
  - [ggrepel](https://github.com/slowkow/ggrepel)


## Processing, quality check and trimming

Raw `bcl` files were converted to `fastq` using `bcl2fastq`:

```bash
bcl2fastq --sample-sheet SampleSheet.csv --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0
```

where `SampleSheet.csv` is the file containing Illumina indexes for each library. Sequencing quality was explored using `fastqc` and looked acceptable for further analysis. Illumina sequencing adaptors were removed from the `fastq` files using `cutadapt`:

```bash
cutadapt -a AGATCGGAAGAGC -o file_trimmed.fastq.gz file.fastq.gz
```

where `file.fastq.gz` is the input sequencing files and the `_trimmed` counterpart is the output.

Sequences for the widom templates were:

- Forward: ATCGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCGAT

- Reverse: ATCGGATGTATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACAGCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCCTCGGCACCGGGATTCTCGAT


## Alignment

Alignments to the widom templates were generated using `bowtie2`:

```bash
bowtie2 -L 20 -x $ref -U $fq | samtools view -b - | samtools sort -T /tmp/$bname -o $bam - && samtools index $bam
```
where depending on whether the trimmed input `$fq` refers to forward or reverse then `$ref` is a `fasta` file containing the forward or the reverse widom template respectively. `samtools sort` and `samtools index` were used to sort and index the resulting alignment `$bam` file.

The aligner `bwa mem` was also tested however unlike `bowtie2` it was not able to provide alignments in the regions close to the start of the template unfortunately:

```bash
bwa mem -t 4 -k 5 -M $ref $fq | samtools view -b - | samtools sort -T /tmp/$bname -o $bam - && samtools index $bam
```


## Counts

Counting of alignment start positions of matched reads to their corresponding templates was performed as follows:

```bash
for bam in *.bam
do
bname=${bam%.bam}
echo $bname
samtools view $bam -F 4 -F 16 | awk '$6 !~ /[IDSH]/' | awk '$10 !~ /[N]/' | awk '{ print $4 + length($10) - 1 }' | sort -k1,1n | uniq -c | awk ' { t = $1; $1 = $2; $2 = t; print; } ' > $bname.txt
done
```

The resulting `txt` files were organised into forward and reverse tables using the following script:

```python
import os

libs = sorted(os.listdir("/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/data/20170627/counts_bowtie2"))

#reverse
seq_rev = "TAGCCTACATATATAGACTGTGCACGGACCTCTGATCCCTCATTAGGGGAACCGCCAATTTTGCGCCCCCTGTCGCGCATGCACGCAAATTCGCCACGATCTCGACAGATGCTGGTTAACTCGCCGGAGCCGTGGCCCTAAGAGCTA"

lib_pos_cnt = {}

for lib in libs:
  if ("Rev" in lib) or ("rev" in lib):
    print lib
    revfile = open("/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/data/20170627/counts_bowtie2/%s" % lib, "r")
    revlines = revfile.readlines()
    revfile.close()
    lib_pos_cnt[lib] = []
    for line in revlines:
      fields = line.split()
      start = int(fields[0])-1
      count = int(fields[1])
      lib_pos_cnt[lib].append((start, seq_rev[start], count))


for lib in lib_pos_cnt:
  for i in range(len(seq_rev)):
    if i not in [e[0] for e in lib_pos_cnt[lib]]:
      lib_pos_cnt[lib].append((i, seq_rev[i], 0))


libs_ordered = ["170614_M00886_0170_000000000-B978R_PolstopFreeRev1_S6.txt", "170621_M00886_0173_000000000-B97B9_PolstopFreerev4_S7.txt", "170621_M00886_0173_000000000-B97B9_PolstopCLrev7_S8.txt", "170614_M00886_0170_000000000-B978R_PolstopCrosslinkRev3_S8.txt"]

ofile = open("/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/tables/20170704_table_rev.txt", "w")

header = "pos\tbase\t"+"\t".join(["free_rev_1", "free_rev_4", "cl_rev_7", "cl_rev_3"])+"\n"

ofile.write(header)

for trip in sorted(lib_pos_cnt[libs_ordered[0]]):
  line="%s\t%s" % (str(trip[0]+1), trip[1])
  for lib in libs_ordered:
    for pos in lib_pos_cnt[lib]:
      if pos[0] == trip[0]:
        line=line+"\t%s" % str(pos[2])
  line=line+"\n"
  ofile.write(line)

ofile.close()



#forward
seq_fwd = "TAGCTCTTAGGGCCACGGCTCCGGCGAGTTAACCAGCATCTGTCGAGATCGTGGCGAATTTGCGTGCATGCGCGACAGGGGGCGCAAAATTGGCGGTTCCCCTAATGAGGGATCAGAGGTCCGTGCACAGTCTATATATGTAGGCTA"

lib_pos_cnt = {}

for lib in libs:
  if ("FW" in lib) or ("Fw" in lib) or ("fw" in lib):
    print lib
    fwdfile = open("/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/data/20170627/counts_bowtie2/%s" % lib, "r")
    fwdlines = fwdfile.readlines()
    fwdfile.close()
    lib_pos_cnt[lib] = []
    for line in fwdlines:
      fields = line.split()
      start = int(fields[0])-1
      count = int(fields[1])
      lib_pos_cnt[lib].append((start, seq_fwd[start], count))


for lib in lib_pos_cnt:
  for i in range(len(seq_fwd)):
    if i not in [e[0] for e in lib_pos_cnt[lib]]:
      lib_pos_cnt[lib].append((i, seq_fwd[i], 0))


libs_ordered = ["170614_M00886_0170_000000000-B978R_PolstopFreeFw13_S7.txt", "170614_M00886_0170_000000000-B978R_PolsttopFreeFw14_S5.txt", "170614_M00886_0170_000000000-B978R_PolstopCrosslinkFW11_S4.txt", "170621_M00886_0173_000000000-B97B9_PolstopCLfw6_S5.txt"]

ofile = open("/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/tables/20170704_table_fwd.txt", "w")

header = "pos\tbase\t"+"\t".join(["free_fwd_13", "free_fwd_14", "cl_fwd_11", "cl_fwd_6"])+"\n"

ofile.write(header)

for trip in sorted(lib_pos_cnt[libs_ordered[0]]):
  line="%s\t%s" % (str(trip[0]+1), trip[1])
  for lib in libs_ordered:
    for pos in lib_pos_cnt[lib]:
      if pos[0] == trip[0]:
        line=line+"\t%s" % str(pos[2])
  line=line+"\n"
  ofile.write(line)

ofile.close()
```


## Analysis

We analysed the resulting `20170704_table_rev.txt` and `20170704_table_fwd.txt` tables using the R programming language:

```r
library(edgeR)
library(ggplot2)
library(data.table)
library(ggrepel)

# Enlarge the view width when printing tables
options(width = 300)

############
############
# reverse ##
############
############

# Load data
data_rev <- read.table("/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/tables/20170704_table_rev.txt", header = T)

# Define a data objects that we can manipulate using edgeR
data_rev_n <- data_rev[,3:6]
rownames(data_rev_n) <- paste(data_rev$pos, data_rev$base, sep = "_")

# Omit last line (full alignments)
data_rev_n <- data_rev_n[1:nrow(data_rev_n)-1,]

# Define DGEList object and normalise counts
group <- c("free", "free", "cl", "cl")
y_rev <- DGEList(counts = data_rev_n, group = group)
libSize <- c(free_rev_1 = 1550593, free_rev_4 = 1440416, cl_rev_7 = 2163278, cl_rev_3 = 1471188)
y_rev$samples$lib.size <- libSize
y_rev <- calcNormFactors(y_rev, method = "none")

# Differential stop analysis
# Set up design matrix for GLM
des <- model.matrix(~ 0 + group, data = y_rev$samples)
colnames(des) <- levels(factor(y_rev$samples$group))
des

# Estimate dispersions
y_rev_glm <- estimateDisp(y_rev, des)
sqrt(y_rev_glm$common.disp) # 0.8827292

# Fit negative bionomial GLM
y_rev_fit <- glmFit(y_rev_glm, des)

# Define matrix of contrasts
my.contrasts <- makeContrasts(clvsfree = cl - free, levels = des)

# Carry out Likelihood ratio tests
lrt_rev_clvsfree <- glmLRT(y_rev_fit, contrast = my.contrasts[,"clvsfree"])


#####################################################
# Show top ranked positions and write them to files #
#####################################################

topTags(lrt_rev_clvsfree, n = 20)
#Coefficient:  1*cl -1*free
#         logFC    logCPM        LR       PValue         FDR
#72_G  3.996800  7.566189 15.717055 7.355794e-05 0.004838427
#98_G  4.029108  8.946845 14.844037 1.167767e-04 0.004838427
#33_T  5.084292  3.846373 14.728724 1.241405e-04 0.004838427
#61_T  5.153896  5.680348 14.605039 1.325596e-04 0.004838427
#59_T  3.513588  5.228024 13.750018 2.088188e-04 0.006036576
#100_T 3.611178  9.103414 13.426622 2.480784e-04 0.006036576
#73_T  3.490560  6.498385 12.775049 3.512734e-04 0.007326560
#70_C  3.501894  7.667640 12.417548 4.253183e-04 0.007402959
#41_C  3.428588  4.147473 12.142655 4.928158e-04 0.007402959
#44_T  3.635746  4.142489 12.089550 5.070520e-04 0.007402959
#31_T  4.680440  2.764020 11.871233 5.700765e-04 0.007566469
#122_C 4.390401 10.427565 11.587454 6.639824e-04 0.008078453
#58_A  3.083959  5.917483 10.569589 1.149633e-03 0.012700460
#40_T  3.940026  5.106270 10.463062 1.217852e-03 0.012700460
#113_T 3.955215  7.909012 10.285149 1.341052e-03 0.013052908
#106_C 2.967821  9.331861  9.613591 1.931426e-03 0.017115885
#34_G  4.957826  6.442025  9.556022 1.992945e-03 0.017115885
#60_T  4.048923  5.597935  9.305663 2.284466e-03 0.018529558
#87_A  3.487985  8.345670  9.037298 2.645264e-03 0.020138457
#42_A  3.063851  3.257297  8.930678 2.804202e-03 0.020138457

y_rev$counts["72_G",] # top1
#free_rev_1 free_rev_4   cl_rev_7   cl_rev_3
#        35         31        725        546
cpm(y_rev)["72_G",]
#free_rev_1 free_rev_4   cl_rev_7   cl_rev_3
#  22.57201   21.52156  335.13954  371.12864

y_rev$counts["98_G",] # top2
#free_rev_1 free_rev_4   cl_rev_7   cl_rev_3
#       110         61       2003       1363
cpm(y_rev)["98_G",]
#free_rev_1 free_rev_4   cl_rev_7   cl_rev_3
#  70.94060   42.34888  925.90966  926.46215

y_rev$counts["145_C",] # bottom
#free_rev_1 free_rev_4   cl_rev_7   cl_rev_3
#     66020      56141     117618      87212
cpm(y_rev)["145_C",]
#free_rev_1 free_rev_4   cl_rev_7   cl_rev_3
#  42577.26   38975.55   54370.27   59279.98

write.table(data.table(as.data.frame(topTags(lrt_rev_clvsfree, n = Inf)), keep.rownames=TRUE)[,.(rn, logFC, PValue, FDR)], file = "/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/tables/20170705_clvsfree_rev.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync -arvuP --remove-source-files /lustre/sblab/martin03/repository/20151125_5fC_nucleosome/tables/20170705_clvsfree_rev.txt martin03@sblab-srv001:/nas/sblab_data1/group_folders/martin03/github/sblab-bioinformatics/projects/20151125_5fC_nucleosome/20170130_robyn/tables/")


################
# Volcano plot #
################

gg <- ggplot(data = data.table(as.data.frame(topTags(lrt_rev_clvsfree, n = Inf)), keep.rownames=TRUE), aes(x = logFC, y = -log10(FDR))) +
geom_point(size = 1) +
theme_bw() +
xlab(expression("log"[2]*"FC")) +
ylab(expression("-log"[10]*"FDR")) +
geom_text_repel(data = data.table(as.data.frame(topTags(lrt_rev_clvsfree)), keep.rownames=TRUE), aes(label = rn), size = 4, force = 1, segment.size = 0)

ggsave('/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_rev_volcano.png', width= 12, height= 12, units= 'cm')
system("rsync -arvuP --remove-source-files /lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_rev_volcano.png martin03@sblab-srv001:/nas/sblab_data1/group_folders/martin03/github/sblab-bioinformatics/projects/20151125_5fC_nucleosome/20170130_robyn/figures/")


################
# Volcano plot # no labels
################

gg <- ggplot(data = data.table(as.data.frame(topTags(lrt_rev_clvsfree, n = Inf)), keep.rownames=TRUE), aes(x = logFC, y = -log10(FDR))) +
geom_point(size = 1) +
theme_bw() +
xlab(expression("log"[2]*"FC")) +
ylab(expression("-log"[10]*"FDR"))

ggsave('/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170710_clvsfree_rev_volcano_nolabels.png', width= 12, height= 12, units= 'cm')


#######
# MDS #
#######

png('/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_rev_mds.png')
p <- plotMDS(y_rev, col = c(rep("black", 2), rep("deepskyblue3", 2)), xlab = "Dimension 1", ylab = "Dimension 2")
dev.off()

system("rsync -arvuP --remove-source-files /lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_rev_mds.png martin03@sblab-srv001:/nas/sblab_data1/group_folders/martin03/github/sblab-bioinformatics/projects/20151125_5fC_nucleosome/20170130_robyn/figures/")


#######
# hca #
#######

hc_rev <- hclust(dist(t(cpm(y_rev))), method = "ward.D")

png('/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_rev_hca.png')
plot(as.dendrogram(hc_rev), horiz = TRUE)
dev.off()

system("rsync -arvuP --remove-source-files /lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_rev_hca.png martin03@sblab-srv001:/nas/sblab_data1/group_folders/martin03/github/sblab-bioinformatics/projects/20151125_5fC_nucleosome/20170130_robyn/figures/")


###########
# FC plot #
###########

data_rev_clvsfree_plot <- data.table(as.data.frame(topTags(lrt_rev_clvsfree, n = Inf)), keep.rownames=TRUE)
temp  <- strsplit(data_rev_clvsfree_plot$rn, "[_]")
temp <- matrix(unlist(temp), ncol=2, byrow=TRUE)
data_rev_clvsfree_plot[, position := temp[,1]]
data_rev_clvsfree_plot[, base := temp[,2]]
data_rev_clvsfree_plot <- data_rev_clvsfree_plot[order(as.numeric(position))]

gg <- ggplot(data = data_rev_clvsfree_plot, aes(x = position, y = 2^(logFC), group = 1)) +
geom_line() +
scale_x_discrete(limits = as.character(data_rev_clvsfree_plot$position), labels = as.character(data_rev_clvsfree_plot$base)) +
xlab("") +
ylab("FC") +
theme_classic() +
geom_vline(xintercept = 33.5, linetype = "dashed") +
ggtitle("reverse")

ggsave("/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_rev_fc.png", width = 50/2.54, height = 10/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_rev_fc.png martin03@sblab-srv001:/nas/sblab_data1/group_folders/martin03/github/sblab-bioinformatics/projects/20151125_5fC_nucleosome/20170130_robyn/figures")


############
# FDR plot #
############

gg <- ggplot(data = data_rev_clvsfree_plot, aes(x = position, y = -log10(FDR), group = 1)) +
geom_line() +
scale_x_discrete(limits = as.character(data_rev_clvsfree_plot$position), labels = as.character(data_rev_clvsfree_plot$base)) +
xlab("") +
ylab(expression("-log"[10]*"FDR")) +
theme_classic() +
geom_vline(xintercept = 33.5, linetype = "dashed") +
ggtitle("reverse")

ggsave("/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_rev_fdr.png", width = 50/2.54, height = 10/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_rev_fdr.png martin03@sblab-srv001:/nas/sblab_data1/group_folders/martin03/github/sblab-bioinformatics/projects/20151125_5fC_nucleosome/20170130_robyn/figures")



############
############
# forward ##
############
############

# Load data
data_fwd <- read.table("/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/tables/20170704_table_fwd.txt", header = T)

# Define a data objects that we can manipulate using edgeR
data_fwd_n <- data_fwd[,3:6]
rownames(data_fwd_n) <- paste(data_fwd$pos, data_fwd$base, sep = "_")

# Omit last line (full alignments)
data_fwd_n <- data_fwd_n[1:nrow(data_fwd_n)-1,]

# Define DGEList object and normalise counts
group <- c("free", "free", "cl", "cl")
y_fwd <- DGEList(counts = data_fwd_n, group = group)
libSize <- c(free_fwd_13 = 1302311, free_fwd_14 = 1088471, cl_fwd_11 = 1203032, cl_fwd_6 = 1340667)
y_fwd$samples$lib.size <- libSize
y_fwd <- calcNormFactors(y_fwd, method = "none")

# Differential stop analysis
# Set up design matrix for GLM
des <- model.matrix(~ 0 + group, data = y_fwd$samples)
colnames(des) <- levels(factor(y_fwd$samples$group))
des

# Estimate dispersions
y_fwd_glm <- estimateDisp(y_fwd, des)
sqrt(y_fwd_glm$common.disp) # 0.7627498

# Fit negative bionomial GLM
y_fwd_fit <- glmFit(y_fwd_glm, des)

# Define matrix of contrasts
my.contrasts <- makeContrasts(clvsfree = cl - free, levels = des)

# Carry out Likelihood ratio tests
lrt_fwd_clvsfree <- glmLRT(y_fwd_fit, contrast = my.contrasts[,"clvsfree"])


#####################################################
# Show top ranked positions and write them to files #
#####################################################

topTags(lrt_fwd_clvsfree, n = 20)
#Coefficient:  1*cl -1*free
#        logFC    logCPM       LR       PValue          FDR
#31_A 7.524328  9.121657 39.42577 3.407696e-10 4.975236e-08
#51_G 7.002973 11.775213 31.27110 2.243950e-08 1.638084e-06
#45_G 5.616321  9.265752 27.13580 1.896533e-07 9.229795e-06
#42_G 5.227707  8.587276 22.41257 2.199294e-06 8.027422e-05
#62_G 4.672637  9.630848 21.36547 3.795450e-06 9.623295e-05
#50_C 4.781604  9.024146 20.94296 4.731643e-06 9.623295e-05
#71_C 4.327585 10.671578 20.86563 4.926551e-06 9.623295e-05
#67_C 4.338667 10.558974 20.73546 5.273038e-06 9.623295e-05
#38_A 5.226182  8.028915 20.28466 6.673360e-06 1.082567e-04
#34_C 5.367924  7.713556 19.99126 7.779709e-06 1.088560e-04
#53_G 4.566168  9.234674 19.89031 8.201480e-06 1.088560e-04
#68_A 4.144215 10.436232 19.06077 1.266213e-05 1.540560e-04
#63_C 4.357065  9.413760 18.79656 1.454291e-05 1.633280e-04
#72_G 5.175756 12.308898 18.01706 2.189342e-05 2.283171e-04
#30_T 6.103931  5.656304 17.43995 2.965282e-05 2.886208e-04
#33_C 4.808182  7.687831 16.97083 3.795857e-05 3.463719e-04
#73_C 4.682588 11.926677 16.81488 4.120888e-05 3.539116e-04
#64_G 4.031202  9.326091 16.47637 4.926010e-05 3.995541e-04
#66_G 3.750558 10.397780 16.15245 5.844282e-05 4.490869e-04
#41_T 4.983983  6.951237 15.98172 6.395696e-05 4.668858e-04

y_fwd$counts["31_A",] # top1
#free_fwd_13 free_fwd_14   cl_fwd_11    cl_fwd_6
#          7           7         464        2440
cpm(y_fwd)["31_A",]
#free_fwd_13 free_fwd_14   cl_fwd_11    cl_fwd_6
#    5.37506     6.43104   385.69215  1819.98960

y_fwd$counts["51_G",] # top2
#free_fwd_13 free_fwd_14   cl_fwd_11    cl_fwd_6
#         37          87        4673       13430
cpm(y_fwd)["51_G",]
#free_fwd_13 free_fwd_14   cl_fwd_11    cl_fwd_6
#   28.41103    79.92863  3884.35220 10017.40179

y_fwd$counts["131_T",] # bottom
#free_fwd_13 free_fwd_14   cl_fwd_11    cl_fwd_6
#       3119      167892      118676       94780
cpm(y_fwd)["131_T",]
#free_fwd_13 free_fwd_14   cl_fwd_11    cl_fwd_6
#   2394.973  154245.726   98647.418   70696.153

write.table(data.table(as.data.frame(topTags(lrt_fwd_clvsfree, n = Inf)), keep.rownames=TRUE)[,.(rn, logFC, PValue, FDR)], file = "/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/tables/20170705_clvsfree_fwd.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync -arvuP --remove-source-files /lustre/sblab/martin03/repository/20151125_5fC_nucleosome/tables/20170705_clvsfree_fwd.txt martin03@sblab-srv001:/nas/sblab_data1/group_folders/martin03/github/sblab-bioinformatics/projects/20151125_5fC_nucleosome/20170130_robyn/tables/")


################
# Volcano plot #
################

gg <- ggplot(data = data.table(as.data.frame(topTags(lrt_fwd_clvsfree, n = Inf)), keep.rownames=TRUE), aes(x = logFC, y = -log10(FDR))) +
geom_point(size = 1) +
theme_bw() +
xlab(expression("log"[2]*"FC")) +
ylab(expression("-log"[10]*"FDR")) +
geom_text_repel(data = data.table(as.data.frame(topTags(lrt_fwd_clvsfree)), keep.rownames=TRUE), aes(label = rn), size = 4, force = 1, segment.size = 0)

ggsave('/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_fwd_volcano.png', width= 12, height= 12, units= 'cm')
system("rsync -arvuP --remove-source-files /lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_fwd_volcano.png martin03@sblab-srv001:/nas/sblab_data1/group_folders/martin03/github/sblab-bioinformatics/projects/20151125_5fC_nucleosome/20170130_robyn/figures/")


################
# Volcano plot # no labels
################

gg <- ggplot(data = data.table(as.data.frame(topTags(lrt_fwd_clvsfree, n = Inf)), keep.rownames=TRUE), aes(x = logFC, y = -log10(FDR))) +
geom_point(size = 1) +
theme_bw() +
xlab(expression("log"[2]*"FC")) +
ylab(expression("-log"[10]*"FDR"))

ggsave('/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170710_clvsfree_fwd_volcano_nolabels.png', width= 12, height= 12, units= 'cm')


#######
# MDS #
#######

png('/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_fwd_mds.png')
p <- plotMDS(y_fwd, col = c(rep("black", 2), rep("deepskyblue3", 2)), xlab = "Dimension 1", ylab = "Dimension 2")
dev.off()

system("rsync -arvuP --remove-source-files /lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_fwd_mds.png martin03@sblab-srv001:/nas/sblab_data1/group_folders/martin03/github/sblab-bioinformatics/projects/20151125_5fC_nucleosome/20170130_robyn/figures/")


#######
# hca #
#######

hc_fwd <- hclust(dist(t(cpm(y_fwd))), method = "ward.D")

png('/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_fwd_hca.png')
plot(as.dendrogram(hc_fwd), horiz = TRUE)
dev.off()

system("rsync -arvuP --remove-source-files /lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_fwd_hca.png martin03@sblab-srv001:/nas/sblab_data1/group_folders/martin03/github/sblab-bioinformatics/projects/20151125_5fC_nucleosome/20170130_robyn/figures/")


###########
# FC plot #
###########

data_fwd_clvsfree_plot <- data.table(as.data.frame(topTags(lrt_fwd_clvsfree, n = Inf)), keep.rownames=TRUE)
temp  <- strsplit(data_fwd_clvsfree_plot$rn, "[_]")
temp <- matrix(unlist(temp), ncol=2, byrow=TRUE)
data_fwd_clvsfree_plot[, position := temp[,1]]
data_fwd_clvsfree_plot[, base := temp[,2]]
data_fwd_clvsfree_plot <- data_fwd_clvsfree_plot[order(as.numeric(position))]

gg <- ggplot(data = data_fwd_clvsfree_plot, aes(x = position, y = 2^(logFC), group = 1)) +
geom_line() +
scale_x_discrete(limits = as.character(data_fwd_clvsfree_plot$position), labels = as.character(data_fwd_clvsfree_plot$base)) +
xlab("") +
ylab("FC") +
theme_classic() +
geom_vline(xintercept = 20.5, linetype = "dashed") +
ggtitle("forward")

ggsave("/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_fwd_fc.png", width = 50/2.54, height = 10/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_fwd_fc.png martin03@sblab-srv001:/nas/sblab_data1/group_folders/martin03/github/sblab-bioinformatics/projects/20151125_5fC_nucleosome/20170130_robyn/figures")


############
# FDR plot #
############

gg <- ggplot(data = data_fwd_clvsfree_plot, aes(x = position, y = -log10(FDR), group = 1)) +
geom_line() +
scale_x_discrete(limits = as.character(data_fwd_clvsfree_plot$position), labels = as.character(data_fwd_clvsfree_plot$base)) +
xlab("") +
ylab(expression("-log"[10]*"FDR")) +
theme_classic() +
geom_vline(xintercept = 20.5, linetype = "dashed") +
ggtitle("forward")

ggsave("/lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_fwd_fdr.png", width = 50/2.54, height = 10/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20151125_5fC_nucleosome/figures/20170705_clvsfree_fwd_fdr.png martin03@sblab-srv001:/nas/sblab_data1/group_folders/martin03/github/sblab-bioinformatics/projects/20151125_5fC_nucleosome/20170130_robyn/figures")
```


## Positiong of nucleosome phase based on polymerase stop experiments for 5fC

The programs and scripts used are described [here](polymerase_stop_nucleosome_position_5fC.md)
