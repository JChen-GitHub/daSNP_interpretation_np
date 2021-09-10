#================================
#Conduct MAGMA, H-MAGMA or nMAGMA
#================================
##Perform gene analysis using MAGMA 
## original MAGMA
#file_1: file_1: 1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.txt
your-working-directory/magma_v1.07b/magma --bfile your-working-directory/MAGMA_Annotation_file/g1000_eur_ref_file/g1000_eur \
--pval your-working-directory/PGC_sumstat/brain_ncSNP_gwas_sumstat/1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.txt N=50000 \
--gene-annot your-working-directory/MAGMA_Annotation_file/brain_ncSNP_results/annotate/gene_direct_annotate/1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.snp_to_gene.genes.annot \
--gene-model snp-wise=mean \
--out your-working-directory/MAGMA_Annotation_file/brain_ncSNP_results/gene-analysis/MAGMA/1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.MAGMA.gene_analysis

##H-MAGMA
#file_1
your-working-directory/magma_v1.07b/magma --bfile /share/inspurStorage/home1/yanganyi/Desktop/MAGMA_Annotation_file/g1000_eur_ref_file/g1000_eur \
--pval your-working-directory/PGC_sumstat/brain_ncSNP_gwas_sumstat/1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.txt N=50000 \
--gene-annot your-working-directory/MAGMA_Annotation_file/brain_ncSNP_results/annotate/H-MAGMA_gene_annotate/Fetal_brain_coding.genes.annot \
--gene-model snp-wise=mean \
--out your-working-directory/MAGMA_Annotation_file/brain_ncSNP_results/gene-analysis/H-MAGMA/Fetal_brain/1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.Fetal.H-MAGMA.gene_analysis

##nMAGMA
#file_1
your-working-directory/magma_v1.07b/magma --bfile your-working-directory/MAGMA_Annotation_file/g1000_eur_ref_file/g1000_eur \
--pval your-working-directory/PGC_sumstat/brain_ncSNP_gwas_sumstat/1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.txt N=50000 \
--gene-annot your-working-directory/MAGMA_Annotation_file/brain_ncSNP_results/annotate/nMAGMA_gene_annotate/Lung/1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.Lung.nMAGMA.genes.annot \
--gene-model snp-wise=mean \
--out your-working-directory/MAGMA_Annotation_file/brain_ncSNP_results/gene-analysis/nMAGMA/Lung/1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.Lung.nMAGMA.gene_analysis


