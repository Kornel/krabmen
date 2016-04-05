Krabmen
=======

### Download data

 - To download the necessery TCGA data please run [scripts/r/download-tcga-mrna.R](scripts/r/download-tcga-mrna.R).

 - This script downloads TCGA data into the [download](download) directory.
 
 - To extract headers from the file (for barcode analysis) please run [scripts/first-line.sh](scripts/first-line.sh). These headers can be used to generate barcode stats.

 - Next run [scripts/r/count.R](scripts/r/count.R) to generate statistics of all barcodes in [results/stats](results/stats). These statistics are a summary of barcodes of the TCGA data - how many tumor vs how many normal cells are in each tumor file. The file count.R has a variable `headers.dir` which should be set to either the normalized or raw data header directory.


### Expression level

Run negative binomial test, H0 = mean expression level for tumor and normal cells are the same

 - Run tests in [scripts/r/binomtest](scripts/r/binomtest): [scripts/r/binomtest/run-test.R](scripts/r/binomtest/run-test.R). This generates p-values for each tumor.

 - This results are stored in [results/bionmialTest/partial](results/bionmialTest/partial).

 - To merge these results into one table, please run [scripts/r/binomtest/merge-test-results.R](scripts/r/binomtest/merge-test-results.R).

 - The results are stored in the results directory [results/bionmialTest](results/bionmialTest), in two formats: long and wide. The wide format with tumors in columns are in [results/bionmialTest/full-table-pvalues.csv](results/bionmialTest/full-table-pvalues.csv).


Additionally, a table with statistics like median and mean expression level should be generated for normalized data.

 - To do so, please run [scripts/r/rsem-normalized/generate-tables.R](scripts/r/rsem-normalized/generate-tables.R). The results are stored in [results/rsem-normalized/partial](results/rsem-normalized/partial/).
 
 - Next run [scripts/r/rsem-normalized/merge-tables.R](scripts/r/rsem-normalized/merge-tables.R) to merge all partial results. The result is stored in [results/rsem-normalized/full-table-long.csv](results/rsem-normalized/full-table-long.csv) and [results/rsem-normalized/full-table.csv](results/rsem-normalized/full-table.csv). Long and wide format respectively.

### Heatmaps:

To generate data for heatmaps run:

 - [scripts/r/heatmap/generate-data.R](scripts/r/heatmap/generate-data.R) - generates data based on the results of the binomial test and the stats of the normalized counts.

 - [scripts/r/heatmap/generate-plots.R](scripts/r/heatmap/generate-plots.R) - generates plots based on the data prepared in the previous step. All plots are generated in [scripts/r/heatmap/plots](scripts/r/heatmap/plots)
