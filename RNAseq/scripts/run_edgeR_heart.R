library(edgeR)
library(ggplot2)
library(readr)
library(data.table)
source('intensityFilter.R')
source('TPM.R')

targets<- read.table('design.txt', header= TRUE, sep= ' ', row.names= 1)

Group<- rep(NA, nrow(targets))
for(i in 1:nrow(targets)){
    g<- rep(NA, ncol(targets))
    for(j in 1:ncol(targets)){
        g[j]<- paste(names(targets)[j], targets[i,j], sep= '_')
    }
    Group[i]<- paste(g, collapse= '.')
}
Group<- as.factor(Group)

## Design matrix:
## -----------------
design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)
row.names(design)<- row.names(targets)

stopifnot(rowSums(design) == 1) ## Each sample has one and only one combination of conditions
colSums(design)  # N. of reps per condition.

# ----------------------------------------
# Preparing data
# ----------------------------------------

raw.data<- read.delim("htseq_heart.matrix")
d <- raw.data[, 2:ncol(raw.data)]
rownames(d) <- raw.data[, 1]

d <- DGEList(counts = d, group = Group)
dim(d)
cpm.d <- cpm(d)

head(cpm.d)

filter<- rowSums(cpm.d > 1) >= 3

d <- d[filter, ]
d <- calcNormFactors(d)


d<- estimateGLMCommonDisp(d, design)
d<- estimateGLMTagwiseDisp(d, design)
fit<- glmFit(d, design)

# ---------------------------------------
# Transcripts per Million
# ---------------------------------------
gene_length<- fread('Mus_musculus.NCBIM37.67.gene_length.txt')
setnames(gene_length, names(gene_length), c('gene_id', 'gene_name', 'len'))

stopifnot(gene_length$gene_id == raw.data$feature_id)

tpm.d<- TPM(raw.data[, 2:ncol(raw.data)], gene_length$len)
tpm.d<- cbind(gene_length[, list(gene_id, gene_name)], tpm.d)[filter]
write.table(tpm.d, 'transcripts_per_million.txt', row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)


# ---------------------------------------
# RPKM
# ---------------------------------------
rpkm.d<- rpkm(raw.data[, 2:ncol(raw.data)], gene_length$len)
rpkm.d<- cbind(gene_length[, list(gene_id, gene_name)], rpkm.d)[filter]
write.table(rpkm.d, 'heart_rpkm.txt', row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)


# ---------------------------------------
# Contrasts
# ---------------------------------------
## MEMO: Contrasts are ordered alphanumerically
## E.g. `idh_yes - idh_no` is tested as  `idh_no - idh_yes`
## hence +ve logFC means `idh_no > idh_yes`


contrasts<- makeContrasts(

    KO_vs_WT= tdg_ko - tdg_wt,
    levels= design
)

## Compute contrasts and write out DE tables
pal<- colorRampPalette(c("white", "lightblue", "yellow", "red"), space = "Lab")

for (cntr in colnames(contrasts)){
    lrt<- glmLRT(fit, contrast= contrasts[, cntr])
    detable<- topTags(lrt, n= nrow(d))$table
    detable$zScore<- localZ(detable$logCPM, detable$logFC, nbins= 20)
    detable$gene_id<- row.names(detable)
    detable<- merge(gene_length[, list(gene_id, gene_name)], detable, by= 'gene_id')

    pdf(paste('maplot', cntr, 'edgeR.pdf', sep= '.'))
    par(las= 1, mgp= c(1.75, 0.5, 0))
    smoothScatter(x= detable$logCPM, y= detable$logFC, xlab= 'logCPM', ylab= 'logFC', main= cntr, colramp= pal, col= 'blue')
    points(x= detable$logCPM, y= detable$logFC, col= ifelse(detable$FDR < 0.05, 'grey60', NA), pch= '+', cex= 0.6)
    points(x= detable$logCPM, y= detable$logFC, col= ifelse(abs(detable$zScore) > 2 & detable$FDR < 0.05, 'grey20', NA), pch= '+', cex= 0.6)
    legend('topright', pch= '+', col= c('grey60', 'grey20'), legend= c(paste('N. FDR < 0.05=', sum(detable$FDR < 0.05)), paste('N. FDR < 0.05 & Z > 2=', sum(detable$FDR < 0.05 & abs(detable$zScore) > 2)) ))
    grid(col= 'grey50')
    dev.off()

    write.table(detable, paste('expr_diff_heart', cntr, 'edgeR.txt', sep= '.'),
        row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
}
