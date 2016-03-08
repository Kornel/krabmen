library(knitr)
library(pheatmap)
library(RColorBrewer)
library(reshape2)

plot.heatmap <- function(h, what.to.plot, threshold, threshold.func, filename) {
  
  palette <- rev(brewer.pal(8,"RdYlGn")) # colorRampPalette(c('green', 'yellow', 'red'))(n = 1000)
  
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
