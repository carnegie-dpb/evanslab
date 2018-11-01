##
## compute TPM from the featureCounts output
##

TPM.YX24 = read.table("YX24.B73.featureCounts.txt", header=TRUE)
TPM.YX24$Counts = TPM.YX24$STAR.B73.YX24.Aligned.sortedByCoord.out.bam
TPM.YX24$Chr = NULL
TPM.YX24$Start = NULL
TPM.YX24$End = NULL
TPM.YX24$Strand = NULL
TPM.YX24$STAR.B73.YX24.Aligned.sortedByCoord.out.bam = NULL

TPM.S364 = read.table("S364.B73.featureCounts.txt", header=TRUE)
TPM.S364$Counts = TPM.S364$STAR.B73.S364.Aligned.sortedByCoord.out.bam
TPM.S364$Chr = NULL
TPM.S364$Start = NULL
TPM.S364$End = NULL
TPM.S364$Strand = NULL
TPM.S364$STAR.B73.S364.Aligned.sortedByCoord.out.bam = NULL

TPM.YX24$RPK = TPM.YX24$Counts/(TPM.YX24$Length/1000)
TPM.S364$RPK = TPM.S364$Counts/(TPM.S364$Length/1000)

TPM.YX24.scale = sum(TPM.YX24$RPK)/1e6
TPM.S364.scale = sum(TPM.S364$RPK)/1e6

TPM.YX24$TPM = TPM.YX24$RPK/TPM.YX24.scale
TPM.S364$TPM = TPM.S364$RPK/TPM.S364.scale

