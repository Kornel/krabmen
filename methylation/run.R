library(readr)
library(dplyr)
library(reshape2)
file <- '../data/methylation/BRCA.methylation.27k.450k.txt.gz'

rawdata <- read_delim(file, delim = '\t')

#Only complete cases
data <- rawdata[complete.cases(rawdata),]

#Short colnames
colnames(data) <- substr(colnames(data), 0, 12)

#Subtypes data for clusters
subtypes <- read_delim('../KRAB_and_TCGA_Subtypes/Data v3 2016.06.06/BRCA.csv', delim = ';')
subtypes <- subtypes[, c('sampleID', 'methylation.Clusters')]
subtypes <- subtypes[complete.cases(subtypes),]

#Transpose probes to columns
probes <- data$id
samples <- colnames(data)

#Remove probe names
data <- data[,-1]

#Transpose
tdata <- as.data.frame(t(data))

#Join with subtypes to get methylation.Cluster
tdata$sampleID <- rownames(tdata)
jdata <- tdata %>% inner_join(subtypes, by = c('sampleID' = 'sampleID'))

colnames(jdata) <- c(probes, c('sampleID', 'cluster'))

#Anova

.append.stat <- function(valueCol, discriminator, stat, name, df) {
  stats <- tapply(valueCol, discriminator, stat)
  names(stats) <- paste0(name, names(stats))
  cbind(df, t(as.data.frame(stats)))
}

result <- data.frame()

for (probe in probes) {
  probeCol <- jdata[, probe]
  clusterCol <- as.character(jdata$cluster)
  # Anova
  a <- aov(probeCol ~ clusterCol)
  
  #TukeyHSD
  t.hsd <- as.data.frame(TukeyHSD(a)$`clusterCol`)
  
  # Convert TukeyHSD result to dataframe
  t.hsd$type <- rownames(t.hsd)
  t.hsd$probe <- probe
  t.hds.df <- dcast(t.hsd, probe ~ type, value.var = c('p adj'))
  colnames(t.hds.df)[-1] <- paste0('p-adj-', colnames(t.hds.df)[-1])
  
  # Append stats: mean, SEM, Q1, Q3
  t.hds.df <- .append.stat(probeCol, clusterCol, mean, 'mean-', t.hds.df)
  t.hds.df <- .append.stat(probeCol, clusterCol, function(x) sqrt(var(x) / length(x)), 'SEM-', t.hds.df)
  t.hds.df <- .append.stat(probeCol, clusterCol, function(x) quantile(x, 1/4), 'Q1-', t.hds.df)
  t.hds.df <- .append.stat(probeCol, clusterCol, function(x) quantile(x, 3/4), 'Q3-', t.hds.df)
  t.hds.df <- .append.stat(probeCol, clusterCol, length, 'count-', t.hds.df)

  result <- rbind(result, t.hds.df)
}

write.table(result, file = 'results/BRCA.csv', row.names = F, sep = ',')

# Map 27k probe <-> genes
library(IlluminaHumanMethylation27k.db)
ProbeToSymbol <- IlluminaHumanMethylation27kSYMBOL
mapped_probes <- mappedkeys(ProbeToSymbol)
mapped_probes.df <- as.data.frame(ProbeToSymbol[mapped_probes]) %>%
  dplyr::rename(symbol27k = symbol, probe27k = probe_id)

result27k <- result %>% left_join(mapped_probes.df, by = c('probe' = 'probe27k'))

# Map 450k probe <-> genes
mapping450k <- read.delim('../data/IlluminaHumanMethylation450k-probes-mapping.csv.gz', sep = ',', col.names = c('probe', 'gene', 'multi'))
mapping450k$gene <- as.character(mapping450k$gene)
s2eg <- org.Hs.egSYMBOL2EG
mapping <- as.data.frame(s2eg[mappedkeys(s2eg)])

probe2symbol <- mapping450k %>%
  inner_join(mapping, by = c('gene' = 'gene_id')) %>%
  dplyr::select(probe, symbol) %>%
  dplyr::rename(probe450k = probe, symbol450k = symbol)

result450k <- result27k %>% left_join(probe2symbol, by = c('probe' = 'probe450k'))

write.table(result450k, file = 'results/BRCA-with-probe-mapping.csv', row.names = F, sep = ',')