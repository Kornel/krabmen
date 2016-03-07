library(knitr)
library(pheatmap)
library(RColorBrewer)
library(reshape2)

heatmap.data <- read.csv('../../results/heatmap/heatmap-data.csv')
h <- heatmap.data[, c('log2median', 'log2mean', 'Gene.Id', 'tumor.name')]  

plot.heatmap <- function(which, h, threshold, filename) {
  palette <- colorRampPalette(c('green', 'yellow', 'red'))(n = 1000)
  
  h.wide <- dcast(h, tumor.name ~ Gene.Id, value.var = paste0('log2', which))
  rownames(h.wide) <- h.wide$tumor.name
  h.wide$tumor.name <- NULL
  h.wide <- t(h.wide)
  h.wide <- data.frame(h.wide)
  #h.wide <- h.wide[complete.cases(h.wide),]
  
  f <- h.wide[abs(rowMeans(h.wide, na.rm = T)) >= threshold,]
  
  dist.rows <- dist(f)
  dist.rows[is.na(dist.rows)] <- max(dist.rows, na.rm = T)
  cluster.rows <- hclust(dist.rows)
  
  dist.cols <- dist(t(f))
  dist.cols[is.na(dist.cols)] <- max(dist.cols, na.rm = T)
  cluster.cols <- hclust(dist.cols)
  
  pheatmap(f, 
           show_colnames = T, 
           cluster_cols = cluster.cols, 
           cluster_rows = cluster.rows, 
           fontsize_row = 10,
           filename = filename,
           cellwidth = 10,
           cellheight = 10,
           legend_breaks =  seq(-5, 5, 1),
           color = palette)
}

today <- format(Sys.time(), "%Y-%m-%d")
results.dir <- sprintf('../../results/heatmap/plots/%s', today)

if (!dir.exists(results.dir)) dir.create(results.dir, recursive = T)

for (threshold in seq(0, 1, 0.1)) {
  plot.heatmap('mean', h, threshold, sprintf('%s/%f-heatmap-mean.png', results.dir, threshold))
  plot.heatmap('median', h, threshold, sprintf('%s/%f-heatmap-median.png', results.dir, threshold))
}
