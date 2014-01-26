for (i in seq(0, 15, 3)) {
  dt<-read.table(paste(i, 'h_pos_dists/', i, 'h_pos_dist_summary.txt', sep=''))
  dt<-dt[2:12,]
  mx<-as.matrix(dt)
  class(mx)<-'numeric'
  mx<-t(mx)
  barplot(mx[2:4,], beside=T, names.arg=mx[1,], xlab="Read Size", ylab="Count", col=c('cyan4', 'darkblue', 'darkorchid4'))
  title(paste(i, ' Hour Reading Frame Offsets by Read Size'))
  legend('topright', c('0 offset', '1 offset', '2 offset'), pch=1, col=c('cyan4', 'darkblue', 'darkorchid4'))
  dev.copy(pdf, paste(i, '_hour_pos_dists.pdf', sep=''))
  dev.off()
}

for (size in seq(26, 36)) {
  mx<-matrix(nrow=4, ncol=0)
  for (i in seq(0, 15, 3)) {
    dt<-read.table(paste(i, 'h_pos_dists/', i, 'h_pos_dist_summary.txt', sep=''))
    dt<-dt[2:12,]
    m<-as.matrix(dt)
    class(m)<-'numeric'
    m<-t(m)
    m[1,size-25]<-i
    mx<-cbind(mx, m[,size-25])  
    barplot(mx[2:4,], beside=T, names.arg=mx[1,], xlab="Hour", ylab="Count", col=c('cyan4', 'darkblue', 'darkorchid4'))
  title(paste(size, 'nt Read Reading Frame Offsets by Hour'))
  legend('topright', c('0 offset', '1 offset', '2 offset'), pch=1, col=c('cyan4', 'darkblue', 'darkorchid4'))
  dev.copy(pdf, paste(size, '_nt_pos_dists.pdf', sep=''))
  dev.off()
  }
}
