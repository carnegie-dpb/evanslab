##
## load data from an alleles file
##

W22 = read.table("W22.d10.Q10.txt", header=T)
W22$reftot = W22$reffor + W22$refrev
W22$alttot = W22$altfor + W22$altrev
