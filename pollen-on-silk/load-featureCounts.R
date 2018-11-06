##
## load featureCounts data and compute FPKM and TPM; merge into single dataframe.
##

YX24.counts = read.table("YX24.B73.featureCounts.txt", header=TRUE)
YX24.counts$Counts = YX24.counts$STAR.B73.YX24.Aligned.sortedByCoord.out.bam
YX24.counts$Chr = NULL
YX24.counts$Start = NULL
YX24.counts$End = NULL
YX24.counts$Strand = NULL
YX24.counts$STAR.B73.YX24.Aligned.sortedByCoord.out.bam = NULL

YX24.counts$RPK = YX24.counts$Counts/(YX24.counts$Length/1000)
YX24.counts.scale = sum(YX24.counts$RPK)/1e6
YX24.counts$TPM = YX24.counts$RPK/YX24.counts.scale

YX24.scale = sum(YX24.counts$Counts)/1e6
YX24.counts$FPKM = YX24.counts$Counts/YX24.scale/(YX24.counts$Length/1000)

######################################################################

S364.counts = read.table("S364.B73.featureCounts.txt", header=TRUE)
S364.counts$Counts = S364.counts$STAR.B73.S364.Aligned.sortedByCoord.out.bam
S364.counts$Chr = NULL
S364.counts$Start = NULL
S364.counts$End = NULL
S364.counts$Strand = NULL
S364.counts$STAR.B73.S364.Aligned.sortedByCoord.out.bam = NULL

S364.counts$RPK = S364.counts$Counts/(S364.counts$Length/1000)
S364.counts.scale = sum(S364.counts$RPK)/1e6
S364.counts$TPM = S364.counts$RPK/S364.counts.scale

S364.scale = sum(S364.counts$Counts)/1e6
S364.counts$FPKM = S364.counts$Counts/S364.scale/(S364.counts$Length/1000)

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

both = merge(both, YX24.counts, by="Gene")
both = merge(both, S364.counts, by="Gene", suffixes = c(".YX24", ".S364"))

######################################################################

