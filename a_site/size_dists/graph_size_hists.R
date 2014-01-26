for ((i in seq(0, 15, 3)) {
  dt<-read.table(paste(i, 'h_foot_size_dist.txt', sep=''));
  colnames(dt)<-c('Size', 'Read.Count');
  barplot(dt$Read.Count, names.arg=dt$Size, xlab='Read Size (nt)', ylab='Count');
  title(paste(i, 'Hour Read Size distribution'));
  dev.copy(pdf, paste(i, '_hour_read_size_hist.pdf', sep=''));
  dev.off();
}
