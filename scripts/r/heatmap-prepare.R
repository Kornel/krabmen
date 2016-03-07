prepare.heatmap.data <- function(gene.table) {
  
  t <- gene.table
  
  t.tumor <- t[ , grep("tumor", colnames(t))]
  t.healthy <- t[ , grep("healthy", colnames(t))]
  
  t.tumor$median = apply(t.tumor, 1, median)
  t.tumor$mean = apply(t.tumor, 1, mean)
  
  t.healthy$median = apply(t.healthy, 1, median)
  t.healthy$mean = apply(t.healthy, 1, mean)
  
  t.tumor <- t.tumor[,c('median', 'mean')]
  t.healthy <- t.healthy[,c('median', 'mean')]
  
  colnames(t.tumor) <- paste0('tumor.', colnames(t.tumor))
  colnames(t.healthy) <- paste0('healthy.', colnames(t.healthy))
  
  t.tumor$Gene.Id <- rownames(t.tumor)
  t.healthy$Gene.Id <- rownames(t.healthy)
  
  t.joined <- t.tumor %>% inner_join(t.healthy, by = c('Gene.Id' = 'Gene.Id'))
  t.joined$log2median <- log2(t.joined$tumor.median/t.joined$healthy.median)
  t.joined$log2mean <- log2(t.joined$tumor.mean/t.joined$healthy.mean)
  
  t.joined$tumor.name <- tumor.name
  
  t.joined <- t.joined[!is.na(t.joined$log2median) & !is.infinite(t.joined$log2median), ]
}