source('subtypes/utils.R')

# PRAD
t.f <- '../../download/mrna-rsem-normalized-2015-11-01/gdac.broadinstitute.org_PRAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2015110100.0.0/PRAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt'
data <- load.tumor(t.f)

s.f <- '../../data/PRAD.clinical.tsv'
o.f <- '../../subtypes/result-PRAD.clinical.csv'
generate.pvalues(data, s.f, o.f)

s.f <- '../../data/Table_S1_PRAD.csv'
o.f <- '../../subtypes/result-TABLE_S1_PRAD.csv'
generate.pvalues(data, s.f, o.f)

# THCA

data <- load.tumor('../../download/mrna-rsem-normalized-2015-11-01/gdac.broadinstitute.org_THCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2015110100.0.0/THCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt')

s.f <- '../../data/THCA.clinical.tsv'
o.f <- '../../subtypes/result-THCA.clinical.csv'
generate.pvalues(data, s.f, o.f)

s.f <- '../../data/THCA_RNA_Subtypes.tsv'
o.f <- '../../subtypes/result-THCA_RNA_Subtypes.csv'
generate.pvalues(data, s.f, o.f)

# BRCA

data <- load.tumor('../../download/mrna-rsem-normalized-2015-11-01/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2015110100.0.0/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt')
s.f <- '../../data/subtypes_BRCA.tsv'
o.f <- '../../subtypes/result-subtypes_BRCA.csv'
generate.pvalues(data, s.f, o.f)

# BLCA

data <- load.tumor('../../download/mrna-rsem-normalized-2015-11-01/gdac.broadinstitute.org_BLCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2015110100.0.0/BLCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt')
s.f <- '../../data/BLCA.clinical.tsv'
o.f <- '../../subtypes/result-BLCA.clinical.csv'
generate.pvalues(data, s.f, o.f)

s.f <- '../../data/BLCA_mRNA_Clusters.tsv'
o.f <- '../../subtypes/result-BLCA_mRNA_Clusters.csv'
generate.pvalues(data, s.f, o.f)

# LUSC

data <- load.tumor('../../download/mrna-rsem-normalized-2015-11-01/gdac.broadinstitute.org_LUSC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2015110100.0.0/LUSC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt')
s.f <- '../../data/LUSC.clinical.basic.tsv'
o.f <- '../../subtypes/result-LUSC.clinical.basic.csv'
generate.pvalues(data, s.f, o.f)

