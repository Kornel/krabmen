Krabmen
=======

### Download data

 - To download the necessery TCGA data please run [scripts/r/download-tcga-mrna.R](scripts/r/download-tcga-mrna.R).

 - This script downloads TCGA data into the [download](download) directory.

 - Next run [scripts/r/count.R](scripts/r/count.R) to generate statistics of the data in [results/stats](results/stats). These statistics are a summary of barcodes of the TCGA data - how many tumor vs how many normal cells are in each tumor file.


### Expression level

Run negative binomial test, H0 = mean expression level for tumor and normal cells are the same

 - Run tests in [scripts/r/binomtest](scripts/r/binomtest): [scripts/r/binomtest/run-test.R](scripts/r/binomtest/run-test.R). This generates p-values for each tumor.

 - This results are stored in [results/bionmialTest/partial](results/bionmialTest/partial).

 - To merge these results into one table, please run [scripts/r/binomtest/merge-test-results.R](scripts/r/binomtest/merge-test-results.R).

 - The results are stored in the results directory [results/bionmialTest](results/bionmialTest), in two formats: long and wide. The wide format with tumors in columns are in [results/bionmialTest/full-table-pvalues.csv](results/bionmialTest/full-table-pvalues.csv).


### Heatmaps:

To generate data for heatmaps run:

 - [scripts/r/heatmap/generate-data.R](scripts/r/heatmap/generate-data.R) - generates data based on the results of the binomial test

 - [scripts/r/heatmap/generate-plots.R](scripts/r/heatmap/generate-plots.R) - generates plots based on the data prepared in the previous step. All plots are generated in [scripts/r/heatmap/plots](scripts/r/heatmap/plots)
