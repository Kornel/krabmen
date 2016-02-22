require(data.table)

h.orig <- read.csv('../../results/heatmap/heatmap-data.csv')

h.orig$fold.change <- h.orig$tumor.median / h.orig$healthy.median
h.orig$log2.fold.change <- log2(h.orig$fold.change)

h <- dcast(setDT(h.orig), Gene.Id ~ tumor.name, value.var = c('tumor.median', 'healthy.median', 'fold.change', 'log2.fold.change'))

cancer.names <- gsub('(.*)_(.*)', '\\2', colnames(h)[-1])

order.cols <- unlist(lapply(x[-1], function(a) {
  if (grepl('healthy', a) > 0) return(1)
  if (grepl('tumor', a) > 0) return(2)
  if (grepl('log2.fold.change', a) > 0) return(4)
  if (grepl('fold.change', a) > 0) return(3)
}))

col.order <- c(1, order(cancer.names) + 1) # First is Gene.Id

h <- (as.data.frame(h)[, col.order])

write.csv(h, '../../results/full-table-median.csv', row.names = F)

# Ugly copy&paste

h.orig <- read.csv('../../results/heatmap/heatmap-data.csv')

h.orig$fold.change <- h.orig$tumor.mean / h.orig$healthy.mean
h.orig$log2.fold.change <- log2(h.orig$fold.change)

h <- dcast(setDT(h.orig), Gene.Id ~ tumor.name, value.var = c('tumor.mean', 'healthy.mean', 'fold.change', 'log2.fold.change'))

cancer.names <- gsub('(.*)_(.*)', '\\2', colnames(h)[-1])

order.cols <- unlist(lapply(x[-1], function(a) {
  if (grepl('healthy', a) > 0) return(1)
  if (grepl('tumor', a) > 0) return(2)
  if (grepl('log2.fold.change', a) > 0) return(4)
  if (grepl('fold.change', a) > 0) return(3)
}))

col.order <- c(1, order(cancer.names) + 1) # First is Gene.Id

h <- (as.data.frame(h)[, col.order])

write.csv(h, '../../results/full-table-mean.csv', row.names = F)
