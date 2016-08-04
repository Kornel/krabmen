# Subtypes

## How to reproduce results

 - Open the project file [/scripts/r/r.Rproj](/scripts/r/r.Rproj) in RStudio.
 - Run the file [/scripts/r/subtypes/generate-p-values.R](/scripts/r/subtypes/generate-p-values.R) 

This script runs for each file in [/KRAB_and_TCGA_Subtypes/Data v3 2016.06.06/](/KRAB_and_TCGA_Subtypes/Data%20v3%202016.06.06/) 
an analysis on the files downloaded from TCGA [/download/mrna-rsem-normalized-2015-11-01/](/download/mrna-rsem-normalized-2015-11-01/).

To download the files proceed as in [/README.md](/README.md).

## About

The script performs an ANOVA and then Tukey HSD test for each subtype for each tumor and stores the results for each tumor and subtype in a separate file. The name of the file is the name of the subtype and it is located in a directory named like the tumor it belongs to.


