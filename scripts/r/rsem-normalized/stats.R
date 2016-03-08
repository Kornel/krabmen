source('./utils.R')

get.full.stats <- function(counts) {
  stats <- counts
  colnames(stats) <- gsub('(.*)\\..*', '\\1', colnames(counts))
  
  stats.healthy <- stats[,grep('healthy', colnames(stats))]
  stats.tumor <- stats[,grep('tumor', colnames(stats))]

  stats$meanHealthy <- apply(stats.healthy, 1, mean)
  stats$meanTumor <- apply(stats.tumor, 1, mean)
  
  stats$medianHealthy <- apply(stats.healthy, 1, median)
  stats$medianTumor <- apply(stats.tumor, 1, median)
  
  stats$Q1healthy <- apply(stats.healthy, 1, function(x) quantile(x, 1/4))
  stats$Q3healthy <- apply(stats.healthy, 1, function(x) quantile(x, 3/4))
  
  stats$Q1tumor <- apply(stats.tumor, 1, function(x) quantile(x, 1/4))
  stats$Q3tumor <- apply(stats.tumor, 1, function(x) quantile(x, 3/4))
  
  stats$SEMhealthy <- apply(stats.healthy, 1, se)
  stats$SEMtumor <- apply(stats.tumor, 1, se)
  
  stats$foldChangeMean <- stats$meanTumor / stats$meanHealthy
  stats$foldChangeMedian <- stats$medianTumor / stats$medianHealthy
  
  stats$log2FoldChangeMean <- log2(stats$foldChangeMean)
  stats$log2FoldChangeMedian <- log2(stats$foldChangeMedian)
    
  stats$Gene.ID <- rownames(stats)
  
  drop.data(stats)
}

