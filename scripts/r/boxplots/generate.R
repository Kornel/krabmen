# WIP

library(ggplot2)
library(readr)
library(ggthemes)
library(dplyr)

asNumeric <- function(x) as.numeric(as.character(x))
factorsNumeric <- function(d) modifyList(d, lapply(d[, sapply(d, is.factor)], asNumeric))

results.dir <- '../../results/boxplot'

genes <- read.csv('boxplots/gene-master-list.csv', sep = '\t')
genes$Gene.ID <- as.character(genes$Gene.ID)

patient.freqs.path <- '../../results/stats/stats-rsem-normalized-1-vs-11.csv'
patient.freqs <- read.csv(patient.freqs.path) %>% filter(normal >= 9)

files.path <- '../../download/mrna-rsem-normalized-2015-11-01'

files <- list.files(path = files.path, pattern = '*-clean.txt', recursive = T)

filtered.files <- c()
  
for (tumor in patient.freqs$tumor.name) {
  contains <- lapply(files, function(x) grepl(tumor, x))
  file <- files[contains == T]
  filtered.files <- c(file, filtered.files)
}

missing.genes <- data.frame(tumor = character(0), gene.id = character(0))

for (file in filtered.files) {
  
  tumor.name <- sub('gdac.broadinstitute.org_(\\w*)\\..*', '\\1', file)
  
  print(sprintf("Reading %s", file))

  data <- read_delim(paste0(files.path, '/', file), delim = '\t')
  
  data$HybRefShort <- sub('\\|.*', '', data$`Hybridization REF`)
  
  joined <- genes %>% left_join(data, by = c('Gene.ID' = 'HybRefShort'))
    
  missing <- (joined %>% filter(is.na(`Hybridization REF`)))
  missing$tumor <- tumor.name

  missing.genes <- rbind(missing.genes, missing[, c('Gene.ID', 'tumor')])

  filtered <- joined %>% filter(!is.na(`Hybridization REF`))
  filtered$Ensembl.ID <- NULL
  filtered$`Hybridization REF` <- NULL

  t <- as.data.frame(t(filtered))
  colnames(t) <- unlist(as.list(t[1, ]))
  t <- t[-1,]
  t <- factorsNumeric(t)
  
  t$type <- factor(unlist(lapply(rownames(t), function(patient.code) {
    elems <- unlist(strsplit(patient.code, '-'))
    sample.vial <- elems[4]
    sample <- as.numeric(substr(sample.vial, 1, 2))
    if (sample == 1) {
      return('tumor')
    } else if (sample == 11) {
      return('healthy')
    } else {
      return('other')
    }
  })))
  
  t <- t %>% filter(type != 'other')
  
  p <- ggplot(t, aes(y = ZNF195, x = type, fill = type)) + geom_boxplot(alpha = 0.5)  +
   scale_fill_manual(values = tableau_color_pal()(3)) + ggtitle(tumor.name)
  
  ggsave(sprintf('%s/%s.png', results.dir, tumor.name), p, dpi = 180)
}

write.csv(missing.genes, file = sprintf('%s/missing-genes.txt', results.dir, tumor.name), row.names = F)
