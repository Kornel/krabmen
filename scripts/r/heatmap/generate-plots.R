library(pheatmap)
library(RColorBrewer)
library(reshape2)

results.dir <- '../../results/heatmap/plots'
heatmap.data <- read.csv(sprintf('%s/heatmap-data.csv', results.dir))

h <- heatmap.data[, c('log2median', 'log2mean', 'Gene.Id', 'tumor.name')]  

plot.heatmap <- function(which, h, threshold, filename) {
  palette <- colorRampPalette(c('green', 'black', 'red'))(n = 1000)
  
  h.wide <- dcast(h, tumor.name ~ Gene.Id, value.var = paste0('log2', which))
  rownames(h.wide) <- h.wide$tumor.name
  h.wide$tumor.name <- NULL
  h.wide <- t(h.wide)
  h.wide <- data.frame(h.wide)
  h.wide <- h.wide[complete.cases(h.wide),]
  pheatmap(h.wide[abs(rowMeans(h.wide)) >= threshold, ], 
           show_colnames = T, 
           cluster_cols = T, 
           cluster_rows = T, 
           fontsize_row = 10,
           filename = filename,
           cellwidth = 10,
           cellheight = 10,
           legend_breaks =  seq(-5, 5, 1),
           color = palette)
}

for (threshold in seq(0, 1, 0.1)) {
  plot.heatmap('mean', h, threshold, sprintf('%s/%f-heatmap-mean.png', results.dir, threshold))
  plot.heatmap('median', h, threshold, sprintf('%s/%f-heatmap-median.png', results.dir, threshold))
}
