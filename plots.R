pdf("estimateh.pdf")

bin<-hexbin(est.raw, est.cmp, xbins=100)
plot(bin, main="Estimate CPM")


bin<-hexbin(est.raw, est.sf, xbins=100)
plot(bin, main="Estimate SF")


bin<-hexbin(est.raw, est.uq, xbins=100)
plot(bin, main="Estimate UQ")


bin<-hexbin(est.raw, est.dms, xbins=100)
plot(bin, main="Estimate DSM")


bin<-hexbin(est.raw, est.udms, xbins=100)
plot(bin, main="Estimate UDMS")

dev.off()

pdf("pvalsh.pdf")

bin<-hexbin(pvals.raw,pvals.cmp, xbins=100)
plot(bin, main="P.values CPM")


bin<-hexbin(pvals.raw, pvals.sf, xbins=100)
plot(bin, main="P.values SF")


bin<-hexbin(pvals.raw, pvals.uq, xbins=100)
plot(bin, main="P.values UQ")


bin<-hexbin(pvals.raw, pvals.dms, xbins=100)
plot(bin, main="P.values DSM")


bin<-hexbin(pvals.raw, pvals.udms, xbins=100)
plot(bin, main="P.values UDMS")

dev.off()
