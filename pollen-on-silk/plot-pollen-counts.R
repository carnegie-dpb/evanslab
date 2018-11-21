################################################################################

## xmin = 1
## xmax = 1e6

## ymin = 1e-3
## ymax = 2e+3

## plot(W22$Var1AT, W22$ratio,
##      xlim=c(xmin,xmax), ylim=c(ymin,ymax), 
##      log="xy", cex=0.5,
##      xlab="YX24 ALT Count", ylab="S364 ALT Count / YX24 ALT Count"
##      )

## lines(c(xmin,xmax), c(W22.ratioBot25,W22.ratioBot25), col="blue", lty=2)
## lines(c(xmin,xmax), c(W22.ratioMedian,W22.ratioMedian), col="blue")
## lines(c(xmin,xmax), c(W22.ratioTop25,W22.ratioTop25), col="blue", lty=2)

## text(W22$Var1AT, W22$ratio, paste(W22$Contig,":",W22$Pos,sep=""), cex=0.6, pos=4, offset=0.2)

################################################################################

xmin = 2
xmax = 2e5

ymin = 2
ymax = 2e5

plot(W22$Var1AT, W22$Var2AT,
     xlim=c(xmin,xmax), ylim=c(ymin,ymax), 
     log="xy", cex=0.3,
     xlab="YX24 ALT Count", ylab="S364 ALT Count",
     main="ALT Reads on W22 per SNP"
     )

## text(W22$Var1AT, W22$Var2AT, paste(W22$Contig,":",W22$Pos,sep=""), cex=0.6, pos=4, offset=0.2)
## text(W22$Var1AT, W22$Var2AT, W22$Gene, cex=0.6, pos=4, offset=0.2)

lines(c(10,xmax), W22.sumRatio*c(10,xmax), col="blue", lty=2, lwd=2)
legend(x="bottomright", inset=0.01,
       c(paste("Ratio of sums =",round(W22.sumRatio,2)),
         "10-folds"),
       lty=c(2,3),
       text.col=c("blue","darkgreen"),
       col=c("blue","darkgreen"),
       bty="n"
       )

lines(c(1,xmax),1e-2*W22.sumRatio*c(1,xmax),col="darkgreen", lty=3)
lines(c(1,xmax),1e-1*W22.sumRatio*c(1,xmax),col="darkgreen", lty=3)
lines(c(1,xmax),1e+1*W22.sumRatio*c(1,xmax),col="darkgreen", lty=3)
lines(c(1,xmax),1e+2*W22.sumRatio*c(1,xmax),col="darkgreen", lty=3)
lines(c(1,xmax),1e+3*W22.sumRatio*c(1,xmax),col="darkgreen", lty=3)
lines(c(1,xmax),1e+4*W22.sumRatio*c(1,xmax),col="darkgreen", lty=3)

## lines(c(1,xmax),W22.ratioBot1*c(1,xmax),col="blue", lty=4)
## lines(c(1,xmax),W22.ratioBot10*c(1,xmax),col="blue", lty=3)
## lines(c(1,xmax),W22.ratioBot25*c(1,xmax),col="blue", lty=2)
## lines(c(1,xmax),W22.ratioMedian*c(1,xmax),col="blue")
## lines(c(1,xmax),W22.ratioTop25*c(1,xmax),col="blue", lty=2)
## lines(c(1,xmax),W22.ratioTop10*c(1,xmax),col="blue", lty=3)
## lines(c(1,xmax),W22.ratioTop1*c(1,xmax),col="blue", lty=4)
## legend(x="bottomright", inset=0.01,
##        c(paste("99 pctl ratio = ",round(W22.ratioTop1,2)),
##          paste("90 pctl ratio = ",round(W22.ratioTop10,2)),
##          paste("75 pctl ratio = ",round(W22.ratioTop25,2)),
##          paste("median ratio =",round(W22.ratioMedian,2)),
##          paste("25 pctl ratio = ",round(W22.ratioBot25,2)),
##          paste("10 pctl ratio = ",round(W22.ratioBot10,2)),
##          paste("1 pctl ratio = ",round(W22.ratioBot1,2))
##          ),
##        text.col=c("blue"),
##        col=c("blue"),
##        lty=c(4,3,2,1,2,3,4),
##        bty="n")

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

## tpmRatio = sum(both$TPM.S364*both$AltFrac.S364)/sum(both$TPM.YX24)

## both.positive = both$Counts.YX24>0 & both$Counts.S364>0

## both.W22.ratioMedian = median(both$TPM.S364[both.positive]/both$TPM.YX24[both.positive])
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
##      main=paste("Gene TPM on B73 for loci with S364 pollen\nmedian ratio = ",round(both.W22.ratioMedian,2),"(blue line)")
##      )

## lines(c(xmin,xmax), both.W22.ratioMedian*c(xmin,xmax), col="blue")

## lines(c(xmin,xmax), 1e-1*both.W22.ratioMedian*c(xmin,xmax), col="red")
## lines(c(xmin,xmax), 1e1*both.W22.ratioMedian*c(xmin,xmax), col="darkgreen")
## lines(c(xmin,xmax), 1e2*both.W22.ratioMedian*c(xmin,xmax), col="darkgreen")
## lines(c(xmin,xmax), 1e3*both.W22.ratioMedian*c(xmin,xmax), col="darkgreen")

## text(xmax, xmax*1e-1*both.W22.ratioMedian, "x0.1", col="red", pos=3)
## text(ymax/both.W22.ratioMedian, ymax, "median", col="blue")
## text(ymax/(1e1*both.W22.ratioMedian), ymax, "x10", col="darkgreen")
## text(ymax/(1e2*both.W22.ratioMedian), ymax, "x100", col="darkgreen")
## text(ymax/(1e3*both.W22.ratioMedian), ymax, "x1000", col="darkgreen")

## ## lines(xmin:xmax-sqrt(xmin:xmax), .10*both.W22.ratioMedian*(xmin:xmax)+sqrt(.10*both.W22.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax+sqrt(xmin:xmax), .10*both.W22.ratioMedian*(xmin:xmax)-sqrt(.10*both.W22.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax-sqrt(xmin:xmax), both.W22.ratioMedian*(xmin:xmax)+sqrt(both.W22.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax+sqrt(xmin:xmax), both.W22.ratioMedian*(xmin:xmax)-sqrt(both.W22.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax-sqrt(xmin:xmax), 10*both.W22.ratioMedian*(xmin:xmax)+sqrt(10*both.W22.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax+sqrt(xmin:xmax), 10*both.W22.ratioMedian*(xmin:xmax)-sqrt(10*both.W22.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax-sqrt(xmin:xmax), 100*both.W22.ratioMedian*(xmin:xmax)+sqrt(100*both.W22.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax+sqrt(xmin:xmax), 100*both.W22.ratioMedian*(xmin:xmax)-sqrt(100*both.W22.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax-sqrt(xmin:xmax), 1000*both.W22.ratioMedian*(xmin:xmax)+sqrt(1000*both.W22.ratioMedian*(xmin:xmax)), col="gray")
## ## lines(xmin:xmax+sqrt(xmin:xmax), 1000*both.W22.ratioMedian*(xmin:xmax)-sqrt(1000*both.W22.ratioMedian*(xmin:xmax)), col="gray")

## ##text(both$TPM.YX24, both$TPM.S364, both$Gene, cex=0.6, pos=4)

################################################################################
