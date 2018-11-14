################################################################################

plot(W22$Var1AF+W22$Var1AR, W22$Var2AF+W22$Var2AR, log="xy", xlab="YX24 ALT Count", ylab="S364 ALT Count", xlim=c(1,2e5), ylim=c(1,5e4), cex=0.5)
lines(c(1,2e5),c(1,5e4),col="blue")

################################################################################

## ymin = 1
## ymax = 1e5

## xmin = 1
## xmax = 1e4

## altLine = 0

## plot(S364$RTnorm+1, S364$ATnorm+1, xlim=c(xmin,xmax), ylim=c(ymin,ymax), cex=0.5,
##      log="xy", 
##      main="S364 reads per B73 gene at YX24 SNP locations on W22",
##      xlab="Avg. REF depth per SNP + 1",
##      ylab="Avg. ALT depth per SNP + 1")

## lines(c(1,xmax), c(1,1), col="black")
## lines(c(1,1), c(1,ymax), col="black")

## xvals = (10^seq(-6,1,1e-3))*xmax/10

## for (f in c(0.001,0.01,0.1,0.5,0.9,0.99,0.999,0.9999)) {
##     yvals = f*xvals/(1-f)
##     lines(xvals+1, yvals+1, col="red", lty=2)
##     text(max(xvals+1), max(yvals+1), f, col="red")
## }
## lines(c(1,1), c(1,ymax), col="red", lty=2)
## text(1, ymax, "1.0", col="red")

## if (altLine>0) {
##     lines(c(xmin,xmax), c(altLine,altLine)+1, col="blue", lty=3)
##     text(xmax, altLine-1, paste("min ALT=",altLine), col="blue", pos=2)
## }

## legend(x="bottomright", inset=0.01, "ALT fraction", text.col="red", lty=2, col="red", bty="n")

## ##text(S364$RTnorm+1, S364$ATnorm+1, S364$Gene, pos=4, cex=0.75, offset=0.2)

################################################################################

## ymin = 10
## ymax = 5e3

## xmin = 4
## xmax = 1e3

## altLine = 0

## countRatio = sum(both$Counts.S364*both$AltFrac.S364)/sum(both$Counts.YX24)

## plot(both$ATnorm.YX24, both$ATnorm.S364,
##      xlim=c(xmin,xmax), ylim=c(ymin,ymax), cex=0.5,
##      log="xy",
##      xlab="YX24 ALT depth per SNP per gene",
##      ylab="S364 ALT depth per SNP per gene",
##      main="ALT depth per B73 gene at YX24 SNP locations on W22"
##      )

## if (altLine>0) {
##     lines(c(altLine,xmax),    c(altLine,altLine), col="blue")
##     lines(c(altLine,altLine), c(altLine,ymax),    col="blue")
##     text(xmax, altLine-1, paste("min ALT=",altLine), col="blue", pos=2)
## }

## lines(c(xmin,xmax), countRatio*c(xmin,xmax), col="blue")
## legend(x="bottomright", bty="n", legend=paste("Adj. count ratio =",round(countRatio,4)), col="blue", text.col="blue", lty=1)

## lines(c(xmin,ymax/(1e1*countRatio)), c(xmin*1e1*countRatio,ymax), col="darkgreen", lty=3)
## lines(c(xmin,ymax/(1e2*countRatio)), c(xmin*1e2*countRatio,ymax), col="darkgreen", lty=3)
## lines(c(xmin,ymax/(1e3*countRatio)), c(xmin*1e3*countRatio,ymax), col="darkgreen", lty=3)

## text(ymax/(1e1*countRatio), ymax, "x10", col="darkgreen")
## text(ymax/(1e2*countRatio), ymax, "x100", col="darkgreen")
## text(ymax/(1e3*countRatio), ymax, "x1000", col="darkgreen")

## text(xmax, xmax*(1e1*countRatio), "x10", col="darkgreen")
## text(xmax, xmax*(1e2*countRatio), "x100", col="darkgreen")
## text(xmax, xmax*(1e3*countRatio), "x1000", col="darkgreen")

## text(both$ATnorm.YX24, both$ATnorm.S364, both$Gene, cex=0.5, pos=4)

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

## tpmRatio = sum(both$TPM.S364*both$AltFrac.S364)/sum(both$TPM.YX24)

## both.positive = both$Counts.YX24>0 & both$Counts.S364>0

## both.ratioMedian = median(both$TPM.S364[both.positive]/both$TPM.YX24[both.positive])
## both.ratioSD = sd(both$TPM.S364[both.positive]/both$TPM.YX24[both.positive])

## both.logRatioMedian = median(log10(both$TPM.S364[both.positive]/both$TPM.YX24[both.positive]))
## both.logRatioSD = sd(log10(both$TPM.S364[both.positive]/both$TPM.YX24[both.positive]))

## ymin = 1e-3
## ymax = 1e4

## xmin = 1e-3
## xmax = 1e5

## plot(both$TPM.YX24, both$TPM.S364*both$AltFrac.S364,
##      xlim=c(xmin,xmax), ylim=c(ymin,ymax), cex=0.5,
##      log="xy",
##      xlab=paste("YX24 TPM (featureCounts)"),
##      ylab=paste("S364 TPM (featureCounts)"),
##      main=paste("Gene TPM on B73 for loci with S364 pollen\nmedian ratio = ",round(both.ratioMedian,2),"(blue line)")
##      )

## lines(c(xmin,xmax), both.ratioMedian*c(xmin,xmax), col="blue")

## lines(c(xmin,xmax), 1e-1*both.ratioMedian*c(xmin,xmax), col="red")
## lines(c(xmin,xmax), 1e1*both.ratioMedian*c(xmin,xmax), col="darkgreen")
## lines(c(xmin,xmax), 1e2*both.ratioMedian*c(xmin,xmax), col="darkgreen")
## lines(c(xmin,xmax), 1e3*both.ratioMedian*c(xmin,xmax), col="darkgreen")

## text(xmax, xmax*1e-1*both.ratioMedian, "x0.1", col="red", pos=3)
## text(ymax/both.ratioMedian, ymax, "median", col="blue")
## text(ymax/(1e1*both.ratioMedian), ymax, "x10", col="darkgreen")
## text(ymax/(1e2*both.ratioMedian), ymax, "x100", col="darkgreen")
## text(ymax/(1e3*both.ratioMedian), ymax, "x1000", col="darkgreen")

## ## lines(xmin:xmax-sqrt(xmin:xmax), .10*both.ratioMedian*(xmin:xmax)+sqrt(.10*both.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax+sqrt(xmin:xmax), .10*both.ratioMedian*(xmin:xmax)-sqrt(.10*both.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax-sqrt(xmin:xmax), both.ratioMedian*(xmin:xmax)+sqrt(both.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax+sqrt(xmin:xmax), both.ratioMedian*(xmin:xmax)-sqrt(both.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax-sqrt(xmin:xmax), 10*both.ratioMedian*(xmin:xmax)+sqrt(10*both.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax+sqrt(xmin:xmax), 10*both.ratioMedian*(xmin:xmax)-sqrt(10*both.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax-sqrt(xmin:xmax), 100*both.ratioMedian*(xmin:xmax)+sqrt(100*both.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax+sqrt(xmin:xmax), 100*both.ratioMedian*(xmin:xmax)-sqrt(100*both.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax-sqrt(xmin:xmax), 1000*both.ratioMedian*(xmin:xmax)+sqrt(1000*both.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax+sqrt(xmin:xmax), 1000*both.ratioMedian*(xmin:xmax)-sqrt(1000*both.ratioMedian*(xmin:xmax)), col="gray")

## ##text(both$TPM.YX24, both$TPM.S364, both$Gene, cex=0.6, pos=4)

################################################################################
