## run all the steps for QTL analysis
step = 10
error.prob = 0.01
n.draws = 64

twinning = est.rf(twinning)
twinning = calc.errorlod(twinning, error.prob=error.prob)
top.errorlod(twinning)
twinning = calc.genoprob(twinning, step=step, error.prob=error.prob)
twinning = sim.geno(twinning, step=step, error.prob=error.prob, n.draws=n.draws)
out = scanone(twinning, verbose=TRUE)
##operm = scanone(twinning, n.perm=100, verbose=TRUE)
##summary(out, perms=operm, alpha=0.05, pvalues=TRUE)

