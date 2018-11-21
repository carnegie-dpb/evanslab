##
## doplot the alt vs ref alleles
##

xmin = 1
xmax = 1e3

ymin = 1e3
ymax = 2e4

plot(genes$SrcRefTot+1, genes$SrcAltTot+1,
     main="B73 pollen on W22 Silk: pollen genes",
     log="xy",
     xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     xlab="W22 SNP REF count + 1", ylab="W22 SNP ALT count + 1",
     pch=1,cex=0.5
     )

text(genes$SrcRefTot+1, genes$SrcAltTot+1, genes$GeneID, pos=4, cex=0.6, offset=0.2)

lines(c(xmin/2,xmax*2), 1*c(xmin/2,xmax*2), col="blue")
lines(c(xmin/2,xmax*2), 1e1*c(xmin/2,xmax*2), col="red")
lines(c(xmin/2,xmax*2), 1e2*c(xmin/2,xmax*2), col="red")
lines(c(xmin/2,xmax*2), 1e3*c(xmin/2,xmax*2), col="red")
lines(c(xmin/2,xmax*2), 1e4*c(xmin/2,xmax*2), col="red")
lines(c(xmin/2,xmax*2), 1e5*c(xmin/2,xmax*2), col="red")

legend("bottomright", cex=0.75, bty="n",
       c("SNPComparer input:",
         "W22.d10.Q10.vcf.gz",
         "remapped_B73_W22.gene.gff3",
         "Zea_mays.B73_RefGen_v4.41.gene.gff3",
         "B73.d10.Q10.vcf.gz",
         paste("sourceAltTotalMin",sourceAltTotalMin),
         paste("sourceAltReadRatioMin",sourceAltReadRatioMin),
         paste("sourceRefFractionMax=0.1",sourceRefFractionMax),
         paste("targetRefFractionMin=0.9",targetRefFractionMin)))

