################################################################################

ratio = (W22$Var2AF+W22$Var2AR)/(W22$Var1AF+W22$Var1AR)

## hist(log2(ratio), breaks=50, main="Histogram of log2(S364/YX24) ALT Counts")

## hist(log2(W22$Var1AF+W22$Var1AR), breaks=20, main="log2(YX24 ALT Counts)", xlim=c(0,20))
hist(log2(W22$Var2AF+W22$Var2AR), breaks=20, main="log2(S364 ALT Counts)", xlim=c(0,20))

################################################################################

## medianLogVal = median(log10(both$Ratio))
## sdLogVal = sd(log10(both$Ratio))

## hist(
##     log10(both$Ratio),
##     breaks=25,
##     main="Ratio of S365/YX24 ALT depth per gene\n(only those without SNPs on B73)",
##     xlab = "log10(S365/YX24)"
## )

## lines(c(medianLogVal,medianLogVal), c(0,par()$yaxp[2]), col="blue")
## lines(c(medianLogVal+sdLogVal,medianLogVal+sdLogVal), c(0,par()$yaxp[2]), col="red")
## lines(c(medianLogVal-sdLogVal,medianLogVal-sdLogVal), c(0,par()$yaxp[2]), col="red")
## lines(c(medianLogVal+2*sdLogVal,medianLogVal+2*sdLogVal), c(0,par()$yaxp[2]), col="red", lty=2)
## lines(c(medianLogVal-2*sdLogVal,medianLogVal-2*sdLogVal), c(0,par()$yaxp[2]), col="red", lty=2)

## legend(x="topright", bty="n",
##        c(paste("median ratio =",round(10^medianLogVal,2)),
##          paste("s.d. ratio =",round(10^sdLogVal,2))),
##        text.col=c("blue","red")
##        )

################################################################################
