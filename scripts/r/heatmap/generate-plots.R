source('./heatmap/plot.R')

heatmap.data <- read.csv('../../results/heatmap/heatmap-data.csv')

today <- format(Sys.time(), "%Y-%m-%d")
results.dir <- sprintf('../../results/heatmap/plots/%s/deseq-normalized', today)

if (!dir.exists(results.dir)) dir.create(results.dir, recursive = T)

for (threshold in seq(0, 1, .1)) {
  
  print(sprintf('Generating heatmaps for threshold %f', threshold))
  
  threshold.func <- function(df) rowMeans(df, na.rm = T)
  mean.path <- sprintf('%s/%f-heatmap-mean-threshold-mean.png', results.dir, threshold)
  median.path <- sprintf('%s/%f-heatmap-median-threshold-mean.png', results.dir, threshold)
  
  plot.heatmap(heatmap.data, 'log2FoldChangeMean', threshold, threshold.func, mean.path)
  plot.heatmap(heatmap.data, 'log2FoldChangeMedian', threshold, threshold.func, median.path)
  
  threshold.func <- function(df) apply(df, 1, median)
  mean.path <- sprintf('%s/%f-heatmap-mean-threshold-median.png', results.dir, threshold)
  median.path <- sprintf('%s/%f-heatmap-median-threshold-median.png', results.dir, threshold)
  
  plot.heatmap(heatmap.data, 'log2FoldChangeMean', threshold, threshold.func, mean.path)
  plot.heatmap(heatmap.data, 'log2FoldChangeMedian', threshold, threshold.func, median.path)
}
