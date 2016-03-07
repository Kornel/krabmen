get.stats <- function(counts) {
  stats <- counts
  colnames(stats) <- gsub('(.*)\\..*', '\\1', colnames(counts))
  
  stats.healthy <- stats[,grep('healthy', colnames(stats))]
  stats.tumor <- stats[,grep('tumor', colnames(stats))]
  
  stats$medianHealthy <- apply(stats.healthy, 1, median)
  stats$medianTumor <- apply(stats.tumor, 1, median)
  
  stats$IQRhealthy <- apply(stats.healthy, 1, IQR)
  stats$IQRtumor <- apply(stats.tumor, 1, IQR)
  
  stats$SEhealthy <- apply(stats.healthy, 1, se)
  stats$SEtumor <- apply(stats.tumor, 1, se)
  
  stats$Gene.ID <- rownames(stats)
  
  subset(stats, select = c('Gene.ID', 'medianHealthy', 'medianTumor', 'IQRhealthy', 'IQRtumor', 'SEhealthy', 'SEtumor'))
}