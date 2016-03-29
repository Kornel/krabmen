library(knitr)
library(pheatmap)
library(RColorBrewer)
library(reshape2)

plot.heatmap <- function(h, what.to.plot, threshold, threshold.func, title, filename) {
  
  h.wide <- dcast(h, tumor.name ~ Gene.Id, value.var = what.to.plot)
  rownames(h.wide) <- h.wide$tumor.name
  h.wide$tumor.name <- NULL
  h.wide <- t(h.wide)
  h.wide <- data.frame(h.wide)
  
  rows <- abs(threshold.func(h.wide)) >= threshold
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
  
  legend.start <- as.integer(min(f, na.rm = T)) - 1
  legend.stop <- as.integer(max(f, na.rm = T)) + 1
  breaks <- seq(legend.start, legend.stop, 1)
  
  zero.pos <- which(breaks == 0)
  greens <- c('#006837', '#1A9850', '#66bd63', '#A6D96A', '#D9EF8B')
  orange <- c('#FEE08B')
  reds <- c('#FDAE61', '#F46D43', '#D73027', '#A50026')
  
  palette.greens <- colorRampPalette(greens, space = 'Lab')(n = zero.pos - 1)
  palette.orange <- orange
  palette.reds <- colorRampPalette(reds, space = 'Lab')(n = length(breaks) - zero.pos + 1)
  palette <- c(palette.greens, palette.orange, palette.reds)
  
  pheatmap(f, 
           show_colnames = T, 
           cluster_cols = cluster.cols, 
           cluster_rows = cluster.rows, 
           fontsize = 10,
           fontsize_row = 10,
           fontsize_col = 10,
           filename = filename,
           cellwidth = 10,
           cellheight = 10,
           treeheight_row = 15,
           treeheight_col = 15,
           color = palette,
           main = title,
           legend_breaks = breaks,
           breaks = breaks)
}

setClass('HeatmapSource', representation(data = 'data.frame', results.dir = 'character'))
setClass('ThresholdMethod', representation(threshold.func = 'function', t.name = 'character'))
