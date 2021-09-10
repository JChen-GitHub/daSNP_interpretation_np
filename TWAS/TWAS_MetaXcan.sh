cd your-working-directory/MetaXcan/software/
conda activate imlabtools

##file_1: 1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.txt
#Cortex
python SPrediXcan.py \
--model_db_path your-working-directory/TWAS/GTEx-V8_elastic_net_eqtls_models/en_Brain_Cortex.db \
--covariance your-working-directory/TWAS/GTEx-V8_elastic_net_eqtls_models/en_Brain_Cortex.txt.gz \
--gwas_folder your-working-directory/PGC_sumstat/brain_ncSNP_gwas_sumstat \
--gwas_file_pattern "1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.txt" \
--snp_column SNP \
--effect_allele_column a1 \
--non_effect_allele_column a2 \
--zscore_column zscore \
--pvalue_column P \
--output_file your-working-directory/TWAS/brain_ncSNP_results/Cortex/1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.Cortex.MetaXcan.csv

