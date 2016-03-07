library(knitr)
library(pheatmap)
library(RColorBrewer)
library(reshape2)

heatmap.data <- read.csv('../../results/heatmap/heatmap-data.csv')

heatmap.data <- do.call(data.frame, lapply(heatmap.data, function(x) replace(x, is.infinite(x), NA)))

plot.heatmap <- function(which, h, threshold, filename) {
  
#   h <- heatmap.data
#   
#   which <- 'Median'
  
  palette <- colorRampPalette(c('green', 'yellow', 'red'))(n = 1000)
  
  h.wide <- dcast(h, tumor.name ~ Gene.Id, value.var = paste0('log2FoldChange', which))
  rownames(h.wide) <- h.wide$tumor.name
  h.wide$tumor.name <- NULL
  h.wide <- t(h.wide)
  h.wide <- data.frame(h.wide)

  rows <- abs(rowMeans(h.wide, na.rm = T)) >= threshold
  rows[is.na(rows)] <- FALSE
  
  f <- h.wide[rows,]
  
  dist.rows <- dist(f)
  dist.rows[is.infinite(dist.rows)] <- NA
  dist.rows[is.na(dist.rows)] <- mean(dist.rows, na.rm = T)
  cluster.rows <- hclust(dist.rows)
  
  dist.cols <- dist(t(f))
  dist.cols[is.infinite(dist.cols)] <- NA
  dist.cols[is.na(dist.cols)] <- mean(dist.cols, na.rm = T)
  cluster.cols <- hclust(dist.cols)
  
  pheatmap(f, 
           show_colnames = T, 
           cluster_cols = cluster.cols, 
           cluster_rows = cluster.rows, 
           fontsize_row = 10,
           filename = filename,
           cellwidth = 10,
           cellheight = 10,
           color = palette)
}

today <- format(Sys.time(), "%Y-%m-%d")
results.dir <- sprintf('../../results/heatmap/plots/%s', today)

if (!dir.exists(results.dir)) dir.create(results.dir, recursive = T)

for (threshold in seq(0, 1, 0.1)) {
  plot.heatmap('Mean', heatmap.data, threshold, sprintf('%s/%f-heatmap-mean.png', results.dir, threshold))
  plot.heatmap('Median', heatmap.data, threshold, sprintf('%s/%f-heatmap-median.png', results.dir, threshold))
}
