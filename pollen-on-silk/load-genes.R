##
## load data from an alleles to genes file
##

## params
sourceAltTotalMin=10
sourceAltReadRatioMin=0.5
sourceRefFractionMax=0.1
targetRefFractionMin=0.9

genes = read.table("pollen-genes-10-0.1-0.9.tsv", header=T, sep="\t")

genes$SrcRefTot = genes$SrcRF + genes$SrcRR
genes$SrcAltTot = genes$SrcAF + genes$SrcAR

genes$GeneID = substring(genes$GeneID,6)
