source("load-pollen-counts.R")

##
## load featureCounts data and compute FPKM and TPM; merge into single dataframe.
##

YX24 = read.table("YX24.B73.featureCounts.txt", header=TRUE)
YX24$Counts = YX24$STAR.B73.YX24.Aligned.sortedByCoord.out.bam
YX24$Chr = NULL
YX24$Start = NULL
YX24$End = NULL
YX24$Strand = NULL
YX24$STAR.B73.YX24.Aligned.sortedByCoord.out.bam = NULL

YX24$RPK = YX24$Counts/(YX24$Length/1000)
YX24.TPMscale = sum(YX24$RPK)/1e6
YX24$TPM = YX24$RPK/YX24.TPMscale

YX24.FPKMscale = sum(YX24$Counts)/1e6
YX24$FPKM = YX24$Counts/YX24.FPKMscale/(YX24$Length/1000)

######################################################################

S364 = read.table("S364.B73.featureCounts.txt", header=TRUE)
S364$Counts = S364$STAR.B73.S364.Aligned.sortedByCoord.out.bam
S364$Chr = NULL
S364$Start = NULL
S364$End = NULL
S364$Strand = NULL
S364$STAR.B73.S364.Aligned.sortedByCoord.out.bam = NULL

S364$RPK = S364$Counts/(S364$Length/1000)
S364.TPMscale = sum(S364$RPK)/1e6
S364$TPM = S364$RPK/S364.TPMscale

S364.FPKMscale = sum(S364$Counts)/1e6
S364$FPKM = S364$Counts/S364.FPKMscale/(S364$Length/1000)

######################################################################

both = merge(YX24, S364, by="Gene", suffixes=c(".YX24",".S364"))

colnames(both)[colnames(both)=="Seq.YX24"] = "Seq"
colnames(both)[colnames(both)=="Start.YX24"] = "Start"
colnames(both)[colnames(both)=="End.YX24"] = "End"
colnames(both)[colnames(both)=="Strand.YX24"] = "Strand"
both$Seq.S364 = NULL
both$Start.S364 = NULL
both$End.S364 = NULL
both$Strand.S364 = NULL

## remove genes with zero counts in both samples
both = both[both$Counts.YX24!=0 | both$Counts.S364!=0,]

## get the combined dataframe of W22 ALT genes, thought to be pollen
pollen = merge(both, YX24, by="Gene")
pollen = merge(pollen, S364, by="Gene", suffixes = c(".YX24", ".S364"))
pollen.ratio = pollen$Counts.S364/pollen$Counts.YX24

######################################################################

