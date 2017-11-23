

# Given a set of genomic intervals, generate an control set with the same CG enrichement profile 

The goal is to generate a set of genomic intervals to use as control, i.e. a set  of the same size as the reference set and that does not overlap with it, which has the same CG enrichment profile. The set of genomic regions will be encoded as a `bed` file. 

If the number of intervals is large, it's computationaly more treatable to split the files in chunks of regions an do the calculations in parallel, merging the control regions at the end. 

In the following, we use the brain TDG KO 5fC regions as example. We use the set of consensus 5fC sites as reference, and the merged (i.e. all detected in the replicates) as sites to avoid when generating the randomized control. We increment the "prohibited" regions
by adding blacklist of unmappable, repetitive, etc. genomic locations.  
In our work, other sets (e.g. heart TDG KO) were generated using the same approach. We make sure that all the regions have the same size by taking the middle of the 5fC region and adding 2000 bp on either side. 

```python
#! /usr/bin/env python
import numpy as np
import pybedtools as pb
from Bio import SeqIO

CHR_SIZES = "PATH_TO_mm9.genome"

bed_5fC_b_ko_merged = pb.BedTool("5fC_merged_KO.bed.gz")
bed_5fC_b_ko = pb.BedTool("5fC_consensus_KO.bed.gz")
LEN_MAX_BLOCK = 1000

def center_regions_bed(bed_file):
    # center of each segment
    long_string = []
    for b in bed_file:
        low = int(0.5 * (b.end + b.start))
        up = int(0.5 * (b.end + b.start))
        string = ("\t").join([str(b.chrom), str(low), str(up)])
        long_string.append(string)
    bed_string = ("\n").join(long_string)
    bed_center = pb.BedTool(bed_string, from_string=True)
    return bed_center

bed_center_5fC_b_ko = center_regions_bed(bed_5fC_b_ko)
sloped_5fC_b_ko = bed_center_5fC_b_ko.slop(l=2000, r=2000, g=CHR_SIZES)
bed_center_5fC_b_ko_merged = center_regions_bed(bed_5fC_b_ko_merged)
sloped_5fC_b_ko_merged = bed_center_5fC_b_ko_merged.slop(l=2000, r=2000, g=CHR_SIZES)

# we are going to split the reference bed file in chuncks
for i in range(0, numb_blocks_b_ko):
    w = i*LEN_MAX_BLOCK
    fname = "5fC_slop_" + str(i) + ".bed"
    _ = pb.BedTool(sloped_5fC_b_ko[w:w+LEN_MAX_BLOCK]).saveas(fname)
w = numb_blocks_b_ko*LEN_MAX_BLOCK
fname = "5fC_slop_" + str(i+1) + ".bed"
_ = pb.BedTool(sloped_5fC_b_ko[w:]).saveas(fname)

# To avoid placing things in "dark" genomic regions, define a blacklist
blacklist = pb.BedTool("mm9-blacklist.bed")

# let's precompute a source of intervals, we do it 6 times to have a larger pool 
# to sample from later on
excluded = sloped_5fC_b_ko_merged.cat(blacklist, postmerge=False)
random_source = sloped_5fC_b_ko.shuffle(g=CHR_SIZES, excl=excluded.fn)
excluded = sloped_5fC_b_ko_merged.cat(random_source.cat(blacklist, postmerge=False), postmerge=False)
REPEAT = 5
for _ in range(REPEAT):
    print(_)
    b = sloped_5fC_b_ko.shuffle(g=CHR_SIZES, excl=excluded.fn)
    xxx = random_source.cat(b, postmerge=False)
    excluded = sloped_5fC_b_ko_merged.cat(xxx.cat(blacklist, postmerge=False), postmerge=False)
    random_source = xxx
    
sloped_rand_cent_5fC_b_ko = random_source 

# we finally split and write the random non-overlaping regions in chunks
# equal in number as the reference set

numb_blocks_b_ko = int(len(sloped_rand_cent_5fC_b_ko)/((REPEAT+1)*LEN_MAX_BLOCK))
for i in range(0, numb_blocks_b_ko):
    w = i*LEN_MAX_BLOCK*(REPEAT+1)
    fname = "rand_non5fC_slop_" + str(i) + ".bed"
    _ = pb.BedTool(sloped_rand_cent_5fC_b_ko[w:w+LEN_MAX_BLOCK*(REPEAT+1)]).saveas(fname)
w = numb_blocks_b_ko*LEN_MAX_BLOCK*(REPEAT+1)
fname = "rand_non5fC_slop_" + str(i+1) + ".bed"
_ = pb.BedTool(sloped_rand_cent_5fC_b_ko[w:]).saveas(fname)
```

## Code to generate the random regions matching the CG enrichment profile or a reference set of regions

The code is [here](https://gitlab.com/guillemportella/fit_randseq_cg_contents)

Run the code in chunks, using the bed files generated above, like this 

```bash
#! /usr/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Give me the max_value for file number, assuming we start at 0 and end at max_value"
    exit
fi

max_value=$1
genome="PATH_TO_genome.fa"

for c in `seq 0 ${max_number}`
do
    bed_ref=ref_beds/5fC_slop_${c}.bed
    bed_rand=rand_beds/rand_non5fC_slop_${c}.bed
    sbatch --wrap "cg_content -i ${genome} -b ${bed_ref} -br ${bed_rand} -max_iter 5000000 -cutoff 0.00010 -bo rand_non5fC_fit_CG_${c}.bed -v"

done 
```

After executions have finished (it takes a long time if you have a large number of reference sites) make sure you reached the cutoff value, then 
concatenate all the rand_fit_CG files into one. 

