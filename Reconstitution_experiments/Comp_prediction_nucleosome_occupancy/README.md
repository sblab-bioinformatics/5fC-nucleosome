
## Computational prediction of nucleosome occupancy

Use van Noort’s sequence-dependent nucleosome occupancy predictor to reproduce the observed in-vitro results.
We have analysed the probability to form nucleosomes in two regions: 200 base pairs around the 5fC site center location, and 200 base pairs flanking the center. How this was done is explained [here](../).

The code used for the calculation can be found [here](nuc_occ_van_noort/). Once compiled, the executable is called `elastic_vnoort`. We used it like this, 

```bash
exe=period_elastic
./$exe -i TDG_KO.Hindbrain.cns4.unique.end.fasta -vnoort  -o occ_vanNoort_end.txt
./$exe -i TDG_KO.Hindbrain.cns4.unique.middle.fasta -vnoort -o occ_vanNoort_middle.txt
./$exe -i TDG_KO.Hindbrain.cns4.unique.start.fasta -vnoort -o occ_vanNoort_start.txt
```


The final figure was prepared using `xmgrace`, I leave [here](xmgrace_file_figure/show_occ_vannoort_regions.agr) the `.agr` file (the format in which `xmgrace` saves its files). 
