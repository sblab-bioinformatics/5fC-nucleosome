# Analysis of Pol-II activity (GRO-seq) in mouse ESC in connection to 5fC

The raw data has been taken from [Wang et al., Nature 2015](https://www.ncbi.nlm.nih.gov/pubmed/26123024)

The GEO for the GRO-seq data is GSE64748. There are 5 datasets (with duplicates) for WT and for TDG-KO.

1.  (WT/KO)-DRB3H --> Time zero of GRO-seq
2.  (WT/KO)-WASH10M --> 10 min after removing the blockage of PolII
3.  (WT/KO)-WASH20M --> 20 min after removing the blockage of PolII
4.  (WT/KO)-WASH30M --> 30 min after removing the blockage of PolII
5.  (WT/KO)-NODRB --> Reference for no blockage of PolII

The protocol is partially based on the original paper, updated based on [Nagari et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5522910/)

### Obtain raw data, perform trimming, alignment and conversion to `BigWig` format.

Described in [this](Groseq_mESC_get_bw.md) document.

### Select transcripts for the metagene analysis, and perform the analysis

Selection of transcripts [here](Groseq_mESC_tx_for_analysis.md)
Generation of metagene profiles [here](Groseq_mESC_metagene_profs.md)

### Plot results

We prepared the figures for publication following this [notebook](groseq_analysis.ipynb) 

