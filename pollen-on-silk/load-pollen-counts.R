##
## load the pollen counts file
## Gene	Seq	Start	End	Strand	Count	RF	RR	AF	AR
##

W22 = read.table(file="W22.SNPFilter.tsv", header=TRUE, sep="\t", fill=TRUE)
## W22 = read.table(file="W22.SNPFilter.purged.tsv", header=TRUE, sep="\t", fill=TRUE)

W22$Gene = substring(W22$Gene,6)
W22$Var1AT = W22$Var1AF + W22$Var1AR
W22$Var2AT = W22$Var2AF + W22$Var2AR
W22$ratio = W22$Var2AT / W22$Var1AT

W22.Var1Perc = quantile(W22$Var1AT, c(.01, .10, .25, .50, .75, .90, .99))
W22.Var2Perc = quantile(W22$Var2AT, c(.01, .10, .25, .50, .75, .90, .99))
W22.ratioPerc = quantile(W22$ratio, c(.01, .10, .25, .50, .75, .90, .99))

W22.ratioBot1 = W22.ratioPerc[1]
W22.ratioBot10 = W22.ratioPerc[2]
W22.ratioBot25 = W22.ratioPerc[3]
W22.ratioMedian= W22.ratioPerc[4]
W22.ratioTop25 = W22.ratioPerc[5]
W22.ratioTop10 = W22.ratioPerc[6]
W22.ratioTop1 = W22.ratioPerc[7]

## these are counted against genes

S364 = read.table(file="S364.W22.SNPReadCounter.0.1.1.1.0.tsv", header=TRUE, sep="\t")
S364$Gene = substring(S364$Gene,6)
S364$RT = S364$RF + S364$RR
S364$AT = S364$AF + S364$AR
S364$AltFrac = S364$AT/(S364$AT + S364$RT)
S364$RTnorm = S364$RT/S364$Count
S364$ATnorm = S364$AT/S364$Count

YX24 = read.table(file="YX24.W22.SNPReadCounter.4.1.9.9.4.tsv", header=TRUE, sep="\t")
YX24$Gene = substring(YX24$Gene,6)
YX24$RT = YX24$RF + YX24$RR
YX24$AT = YX24$AF + YX24$AR
YX24$AltFrac = YX24$AT/(YX24$AT + YX24$RT)
YX24$RTnorm = YX24$RT/YX24$Count
YX24$ATnorm = YX24$AT/YX24$Count
