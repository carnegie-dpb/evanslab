ymin = 1
ymax = 1e5

xmin = 1
xmax = 1e4

altLine = 4
aboveAltLine = both$ATnorm.S364>altLine
medianVal = median(log10(both$ATnorm.S364[aboveAltLine]/both$ATnorm.YX24[aboveAltLine]))
sdVal = sd(log10(both$ATnorm.S364[aboveAltLine]/both$ATnorm.YX24[aboveAltLine]))

################################################################################

plot(S364$RTnorm+1, S364$ATnorm+1, xlim=c(xmin,xmax), ylim=c(ymin,ymax), cex=0.5,
     log="xy", 
     main="S364 reads per B73 gene at YX24 SNP locations on W22\n(only those without SNPs on B73)",
     xlab="Avg. REF depth per SNP + 1",
     ylab="Avg. ALT depth per SNP + 1")

## lines(c(1,xmax), c(1,1), col="black")
## lines(c(1,1), c(1,ymax), col="black")

xvals = (10^seq(-6,1,1e-3))*xmax/10

for (f in c(0.001,0.01,0.1,0.5,0.9,0.99,0.999,0.9999)) {
    yvals = f*xvals/(1-f)
    lines(xvals+1, yvals+1, col="red", lty=2)
    text(max(xvals+1), max(yvals+1), f, col="red")
}
lines(c(1,1), c(1,ymax), col="red", lty=2)
text(1, ymax, "1.0", col="red")

lines(c(xmin,xmax), c(altLine,altLine)+1, col="blue")
text(xmax, altLine-1, paste("ALT=",altLine), col="blue", pos=2)

legend(x="bottomright", inset=0.01, "ALT fraction", text.col="red", lty=2, col="red", bty="n")

##text(S364$RTnorm+1, S364$ATnorm+1, S364$GeneID, pos=4, cex=0.6, offset=0.2)

################################################################################

## plot(both$ATnorm.YX24[aboveAltLine], both$ATnorm.S364[aboveAltLine],
##      xlim=c(4,1e6), ylim=c(4,1e5), cex=1.0,
##      log="xy",
##      xlab="Avg. YX24 ALT depth per SNP per gene",
##      ylab="Avg. S364 ALT depth per SNP per gene",
##      main="ALT reads per B73 gene at YX24 SNP locations on W22\n(only those without SNPs on B73)"
##      )

## lines(c(xmin,xmax), (10^medianVal)*c(xmin,xmax), col="blue")
## lines(c(xmin,xmax), (10^(medianVal+sdVal))*c(xmin,xmax), col="red")
## lines(c(xmin,xmax), (10^(medianVal-sdVal))*c(xmin,xmax), col="red")
## lines(c(xmin,xmax), (10^(medianVal+2*sdVal))*c(xmin,xmax), col="red", lty=2)
## lines(c(xmin,xmax), (10^(medianVal-2*sdVal))*c(xmin,xmax), col="red", lty=2)
## legend(x="bottomright", lty=1, bty="n",
##        c(paste("median ratio=",round(10^medianVal,2)), paste("std.dev.=",round(10^sdVal,2))),
##        col=c("blue","red"),
##        text.col=c("blue","red")
##        )

## text(both$ATnorm.YX24[aboveAltLine], both$ATnorm.S364[aboveAltLine], both$GeneID[aboveAltLine], cex=0.5, pos=4)

################################################################################

## hist(
##     log10(both$ATnorm.S364[aboveAltLine]/both$ATnorm.YX24[aboveAltLine]),
##     breaks=25,
##     main="Ratio of S365/YX24 ALT depth per gene\n(only those without SNPs on B73)",
##     xlab = "log10(S365/YX24)"
## )

## lines(c(medianVal,medianVal), c(0,par()$yaxp[2]), col="blue")
## lines(c(medianVal+sdVal,medianVal+sdVal), c(0,par()$yaxp[2]), col="red")
## lines(c(medianVal-sdVal,medianVal-sdVal), c(0,par()$yaxp[2]), col="red")
## lines(c(medianVal+2*sdVal,medianVal+2*sdVal), c(0,par()$yaxp[2]), col="red", lty=2)
## lines(c(medianVal-2*sdVal,medianVal-2*sdVal), c(0,par()$yaxp[2]), col="red", lty=2)

## legend(x="topright", bty="n",
##        c(paste("median=",round(10^medianVal,4)), paste("s.d.=",round(10^sdVal,4))),
##        text.col=c("blue","red")
##        )

################################################################################
