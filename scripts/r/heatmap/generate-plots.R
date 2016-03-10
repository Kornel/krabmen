source('./heatmap/plot.R')

deseq.heatmap.data <- read.csv('heatmap/data/heatmap-data-deseq-normalized.csv')
rsem.heatmap.data <- read.csv('heatmap/data/heatmap-data-rsem-normalized.csv')

today <- format(Sys.time(), "%Y-%m-%d")

deseq.results.dir <- sprintf('heatmap/plots/%s/deseq-normalized', today)
if (!dir.exists(deseq.results.dir)) dir.create(deseq.results.dir, recursive = T)

rsem.results.dir <- sprintf('heatmap/plots/%s/rsem-normalized', today)
if (!dir.exists(rsem.results.dir)) dir.create(rsem.results.dir, recursive = T)

deseq <- new('HeatmapSource', data = deseq.heatmap.data, results.dir = deseq.results.dir)
rsem <- new('HeatmapSource', data = rsem.heatmap.data, results.dir = rsem.results.dir)
sources <- c(deseq, rsem)

thresh.mean <- new('ThresholdMethod', threshold.func = function(df) rowMeans(df, na.rm = T), t.name = 'mean')
thresh.median <- new('ThresholdMethod', threshold.func = function(df) apply(df, 1, function(x) median(x, na.rm = T)), t.name = 'median')
thresh.methods <- c(thresh.median, thresh.mean)

for (source.data in sources) {
  for (threshold in seq(0, 1, .1)) {
    for (thresh.method in thresh.methods) {
      print(sprintf('Generating heatmaps for threshold %f, method %s, at %s', 
                    threshold, 
                    thresh.method@t.name,
                    source.data@results.dir))
      
      res.dir <- source.data@results.dir
      data <- source.data@data
      thresh.func <- thresh.method@threshold.func
      thresh.name <- thresh.method@t.name
      
      mean.path <- sprintf('%s/%f-heatmap-mean-threshold-%s.png', res.dir, threshold, thresh.name)
      median.path <- sprintf('%s/%f-heatmap-median-threshold-%s.png', res.dir, threshold, thresh.name)
      
      title.mean <- sprintf('log2(mean tumor / mean normal)')
      title.median <- sprintf('log2(median tumor / median normal)')
      
      plot.heatmap(data, 'log2FoldChangeMean', threshold, thresh.func, title.mean, mean.path)
      plot.heatmap(data, 'log2FoldChangeMedian', threshold, thresh.func, title.median, median.path)
    }
  }
}