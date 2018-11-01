##
##
##

both = merge(YX24, S364, by="Gene", suffixes=c(".YX24",".S364"))

colnames(both)[colnames(both)=="Seq.YX24"] = "Seq"
colnames(both)[colnames(both)=="Start.YX24"] = "Start"
colnames(both)[colnames(both)=="End.YX24"] = "End"
colnames(both)[colnames(both)=="Strand.YX24"] = "Strand"
both$Seq.S364 = NULL
both$Start.S364 = NULL
both$End.S364 = NULL
both$Strand.S364 = NULL

both = merge(both, TPM.YX24, by="Gene")
both = merge(both, TPM.S364, by="Gene", suffixes = c(".YX24", ".S364"))
