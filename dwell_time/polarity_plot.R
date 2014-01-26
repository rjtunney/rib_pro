dt=read.csv('plot_data.csv')
dt2<-cbind(dt$Amino.Acid[!duplicated(dt$Amino.Acid)],seq(1:20))
dt2<-dt2[order(dt2[,1]),]
stripchart(dt$Dwell.Time~factor(dt$Amino.Acid), vertical=T, at=dt2[,2], pch=15, col='white')
columns = dt$Amino.Acid[!duplicated(dt$Amino.Acid)]
for (i in 1:20) { 
	aadt = dt[dt$Amino.Acid==columns[i]] 
	points(rep(i,nrow(aadt)), aadt$Dwell.Time, col=ifelse(aadt$Polarity[1]=='Nonpolar', 'blue', 'red'), pch=21, bg='white')
}
pdf(Polarity_Plot.pdf)
