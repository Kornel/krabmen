library(readr)
source('../r/utils.R')

d.raw <- read_delim('../../data/methylation/BRCA.methylation.27k.450k.txt.gz', delim='\t')

rownames(d.raw) <- d.raw$id

# ID is now in the rowname, instead of a regular column
d <- d.raw[, -1]
types <- patient.code.to.type(colnames(d))

# Remove 'other' types
d <- d[, types != "other"]
types <- patient.code.to.type(colnames(d))

# Allocate result dataframe
values <- dim(d)[1]
result <- data.frame(id = character(values), pvalue = numeric(values), stringsAsFactors=FALSE)

for (r in 1:nrow(d)) {
  a <- aov(unlist(d[r,]) ~ types)
  res <- unlist(summary(a))
  id <- rownames(d)[r]
  pval <- res['Pr(>F)1']
  
  result[r, ] <- c(id, pval)
}

# Save results
dir.create('../../scripts/r/methylation/results', showWarnings = F)
write.csv(result, '../../scripts/r/methylation/results/BRCA.csv', row.names = F)