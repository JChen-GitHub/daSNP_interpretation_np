library(limma)
library(stringr)
### gtex expression data process
### gtex data download： https://xenabrowser.net/datapages/?dataset=gtex_RSEM_gene_tpm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
## uniform gene symbol
## expression value cutoff

### gtex data process
gtex_pheno <- read.table('your-working-directory/GTEX/GTEX_phenotype', sep='\t', header=T, row.names = 1 )
rownames(gtex_pheno) <- gsub('-','.', rownames(gtex_pheno) )

grch38 <- read.table('your-working-directory/Homo_sapiens.GRCh38.87.txt', header = T)
rownames(grch38) <- grch38$ensg

gtex_exp <- read.table('your-working-directory/GTEX/gtex_RSEM_gene_tpm',sep='\t',header = T, row.names = 1)

rownames(gtex_exp) <- str_split_fixed(rownames(gtex_exp), '\\.', 2)[,1]

gtex_exp <- gtex_exp[ intersect(rownames(grch38), rownames(gtex_exp)), ]
gtex_exp$gene_symbol <- grch38[rownames(gtex_exp),'symbol']
gtex_exp <- gtex_exp[ !duplicated(gtex_exp$gene_symbol) , ]
rownames(gtex_exp) <- gtex_exp$gene_symbol

gtex_exp <- gtex_exp[, -7863]
gtex_exp <- round(2^gtex_exp,3)-0.001


#  >1 TPM in at least 90% samples--Genes were filtered to include only those with transcripts per million reads (TPM) > 1 in at least 90% of samples,
gtex_brain <- gtex_exp[, intersect(rownames(gtex_pheno)[gtex_pheno$body_site_detail..SMTSD.=='Brain - Frontal Cortex (BA9)'], colnames(gtex_exp)) ]
gtex_brain_exp <- gtex_brain[apply(gtex_brain, 1, function(x){ length(which(x>1)) > length(x)*0.9 }),]

gtex_liver <- gtex_exp[, intersect(rownames(gtex_pheno)[gtex_pheno$body_site_detail..SMTSD.=='Liver'], colnames(gtex_exp)) ]
gtex_Liver_exp <- gtex_liver[apply(gtex_liver, 1, function(x){ length(which(x>1)) > length(x)*0.9 }),]


gtex_lung <- gtex_exp[, intersect(rownames(gtex_pheno)[gtex_pheno$body_site_detail..SMTSD.=='Lung'], colnames(gtex_exp)) ]
gtex_Lung_exp <- gtex_lung[apply(gtex_lung, 1, function(x){ length(which(x>1)) > length(x)*0.9 }),]


gtex_stomach <- gtex_exp[, intersect(rownames(gtex_pheno)[gtex_pheno$body_site_detail..SMTSD.=='Stomach'], colnames(gtex_exp)) ]
gtex_stomach_exp <- gtex_stomach[apply(gtex_stomach, 1, function(x){ length(which(x>1)) > length(x)*0.9 }),]

gtex_intestine <- gtex_exp[, intersect(rownames(gtex_pheno)[gtex_pheno$X_primary_site=='Small Intestine'], colnames(gtex_exp)) ]
gtex_Intestine_exp <- gtex_intestine[apply(gtex_intestine, 1, function(x){ length(which(x>1)) > length(x)*0.9 }),]

gtex_hippocampus <- gtex_exp[, intersect(rownames(gtex_pheno)[gtex_pheno$body_site_detail..SMTSD.%in%c('Brain - Hippocampus')], colnames(gtex_exp)) ]
gtex_hippa_exp <- gtex_hippocampus[apply(gtex_hippocampus, 1, function(x){ length(which(x>1)) > length(x)*0.9 }),]


save(gtex_brain_exp, gtex_Liver_exp, gtex_Lung_exp, gtex_Intestine_exp, gtex_hippa_exp, file='your-working-directory/WGCNA/gtex_wgcna.RData')





