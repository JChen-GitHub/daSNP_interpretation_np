##Annotate SNPs to genes
#file_1: 1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.txt
your-working-directory/magma_v1.07b/magma --annotate \
--snp-loc your-working-directory/MAGMA_Annotation_file/LOC_file/SNPLOC_FILE/brain_ncSNP/1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.SNPLOC.txt \
--gene-loc your-working-directory/MAGMA_Annotation_file/LOC_file/GENELOC_FILE/protein_coding_gene.symbol.loc \
--out your-working-directory/MAGMA_Annotation_file/brain_ncSNP_results/annotate/gene_direct_annotate/1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.snp_to_gene

##Annotate SNPs to REs
#file_1: 1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.txt
your-working-directory/magma_v1.07b/magma --annotate \
--snp-loc your-working-directory/MAGMA_Annotation_file/LOC_file/SNPLOC_FILE/brain_ncSNP/1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.SNPLOC.txt \
--gene-loc your-working-directory/MAGMA_Annotation_file/LOC_file/RELOC_FILE/RELOC_NOGL_hg19.txt \
--out your-working-directory/MAGMA_Annotation_file/brain_ncSNP_results/annotate/re_annotate/1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.snp_to_re

