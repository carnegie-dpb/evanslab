################################################################################

## ymin = 1
## ymax = 1e5

## xmin = 1
## xmax = 1e4

## altLine = 4

## plot(S364$RTnorm+1, S364$ATnorm+1, xlim=c(xmin,xmax), ylim=c(ymin,ymax), cex=0.5,
##      log="xy", 
##      main="S364 reads per B73 gene at YX24 SNP locations on W22\n(only those w/o SNPs on B73)",
##      xlab="Avg. REF depth per SNP + 1",
##      ylab="Avg. ALT depth per SNP + 1")

## ## lines(c(1,xmax), c(1,1), col="black")
## ## lines(c(1,1), c(1,ymax), col="black")

## xvals = (10^seq(-6,1,1e-3))*xmax/10

## for (f in c(0.001,0.01,0.1,0.5,0.9,0.99,0.999,0.9999)) {
##     yvals = f*xvals/(1-f)
##     lines(xvals+1, yvals+1, col="red", lty=2)
##     text(max(xvals+1), max(yvals+1), f, col="red")
## }
## lines(c(1,1), c(1,ymax), col="red", lty=2)
## text(1, ymax, "1.0", col="red")

## if (altLine>0) {
##     lines(c(xmin,xmax), c(altLine,altLine)+1, col="blue")
##     text(xmax, altLine-1, paste("min ALT=",altLine), col="blue", pos=2)
## }

## legend(x="bottomright", inset=0.01, "ALT fraction", text.col="red", lty=2, col="red", bty="n")

## ## text(S364$RTnorm+1, S364$ATnorm+1, S364$Gene, pos=4, cex=0.6, offset=0.2)

################################################################################

## ymin = 1
## ymax = 1e5

## xmin = 1
## xmax = 1e6

## altLine = 4

## medianLogVal = median(log10(both$Ratio))
## sdLogVal = sd(log10(both$Ratio))

## plot(both$ATnorm.YX24, both$ATnorm.S364,
##      xlim=c(xmin,xmax), ylim=c(ymin,ymax), cex=0.5,
##      log="xy",
##      xlab="Avg. YX24 ALT depth per SNP per gene",
##      ylab="Avg. S364 ALT depth per SNP per gene",
##      main="ALT reads per B73 gene at YX24 SNP locations on W22\n(only those without SNPs on B73)"
##      )

## if (altLine>0) {
##     lines(c(altLine,xmax),    c(altLine,altLine), col="blue")
##     lines(c(altLine,altLine), c(altLine,ymax),    col="blue")
##     text(xmax, altLine-1, paste("min ALT=",altLine), col="blue", pos=2)
## }

## lines(c(xmin,xmax), (10^medianLogVal)*c(xmin,xmax), col="red")
## lines(c(xmin,xmax), (10^(medianLogVal+sdLogVal))*c(xmin,xmax), col="red", lty=2)
## lines(c(xmin,xmax), (10^(medianLogVal-sdLogVal))*c(xmin,xmax), col="red", lty=2)

## lines(c(xmin,ymax/1e1), c(ymin*1e1,ymax), col="darkgreen", lty=3)
## lines(c(xmin,ymax/1e2), c(ymin*1e2,ymax), col="darkgreen", lty=3)
## lines(c(xmin,ymax/1e3), c(ymin*1e3,ymax), col="darkgreen", lty=3)

## text(ymax/1e1, ymax, "x10", col="darkgreen")
## text(ymax/1e2, ymax, "x100", col="darkgreen")
## text(ymax/1e3, ymax, "x1000", col="darkgreen")


## legend(x="bottomright",  bty="n",
##        c(paste("median ratio =",round(10^medianLogVal,2)),
##          paste("s.d. ratio =",round(10^sdLogVal,2))),
##        col=c("red","red"),
##        text.col=c("red","red"),
##        lty=c(1,2)
##        )

## ## text(both$ATnorm.YX24, both$ATnorm.S364, both$Gene, cex=0.5, pos=4)

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

ymin = 10
ymax = 1e4

xmin = 1
xmax = 1e3

plot(both$Counts.YX24+1, both$Counts.S364+1,
     xlim=c(xmin,xmax), ylim=c(ymin,ymax), cex=0.5,
     log="xy",
     xlab=paste("YX24 counts+1 on B73"),
     ylab=paste("S364 counts+1 on B73"),
     main="Genes for which reads are (almost) entirely pollen"
     )

lines(c(xmin,xmax), 1e-2*c(xmin,xmax), col="red")
lines(c(xmin,xmax), 1e-1*c(xmin,xmax), col="red")
lines(c(xmin,xmax), 1e0*c(xmin,xmax), col="blue")
lines(c(xmin,xmax), 1e1*c(xmin,xmax), col="darkgreen")
lines(c(xmin,xmax), 1e2*c(xmin,xmax), col="darkgreen")
lines(c(xmin,xmax), 1e3*c(xmin,xmax), col="darkgreen")

text(ymax/1e1, ymax, "x10", col="darkgreen")
text(ymax/1e2, ymax, "x100", col="darkgreen")
text(ymax/1e3, ymax, "x1000", col="darkgreen")


text(both$Counts.YX24+1, both$Counts.S364+1, both$Gene, cex=0.6, pos=4)

################################################################################
