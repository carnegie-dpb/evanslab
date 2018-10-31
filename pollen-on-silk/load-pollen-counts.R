##
## load the pollen counts file
##

S364 = read.table(file="S364.W22.SNPReadCounter.tsv", header=TRUE, sep="\t")

S364$GeneID = substring(S364$GeneID,6)
S364$RT = S364$RF + S364$RR
S364$AT = S364$AF + S364$AR
S364$AltFrac = S364$AT/(S364$AT + S364$RT)
S364$RTnorm = S364$RT/S364$Count
S364$ATnorm = S364$AT/S364$Count


YX24 = read.table(file="YX24.W22.Q10.SNPReadCounter.tsv", header=TRUE, sep="\t")

YX24$GeneID = substring(YX24$GeneID,6)
YX24$RT = YX24$RF + YX24$RR
YX24$AT = YX24$AF + YX24$AR
YX24$AltFrac = YX24$AT/(YX24$AT + YX24$RT)
YX24$RTnorm = YX24$RT/YX24$Count
YX24$ATnorm = YX24$AT/YX24$Count

both = merge(YX24,S364,by="GeneID",suffixes=c(".YX24",".S364"))

