####################################################################

xmin = 1
xmax = 1e6

ymin = 1
ymax = 1e6

plot(both$Counts.YX24, both$Counts.S364,
     xlab="YX24 counts", ylab="S364 counts",
     main="Gene featureCounts on B73",
     xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     cex=0.5, log="xy")

points(pollen$Counts.YX24, pollen$Counts.S364,
       pch=16, cex=0.6, col="darkred")

lines(c(xmin,xmax), 0.10*c(xmin,xmax), col="blue")
lines(c(xmin,xmax), 0.20*c(xmin,xmax), col="blue")
lines(c(xmin,xmax), 0.30*c(xmin,xmax), col="blue")
lines(c(xmin,xmax), median(pollen.Ratio)*c(xmin,xmax), col="red")
legend("bottomright",
       c(paste("median pollen ratio =",round(median(pollen.Ratio),2)),
         "ratio = 0.3",
         "ratio = 0.1"
         ),
       lty=1,
       col=c("red","blue","blue"),
       text.col=c("red","blue","blue"),
       )

## text(both$Counts.YX24, both$Counts.S364, both$Gene,
##      cex = 0.5, pos=4, offset=0.2)

#####################################################################################################################

## hist(log2(pollen.ratio), breaks=50,
##      main=paste("median=",round(median(pollen.ratio),2),"log2(median)=",round(log2(median(pollen.ratio)),2)))

## lines(c(log2(median(pollen.ratio)),log2(median(pollen.ratio))), c(0,1000), col="darkred", lwd=2)

#####################################################################################################################

## hist(log2(pollen$Counts.YX24), breaks=25)

#####################################################################################################################
