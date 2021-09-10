# daSNP_interpretation_np

####################################################################
# I. Construction of TRN (Transcriptional Regulatory Networks)

# Step1 Genome scan of TF motifs.

## s1.1 determine the genome (GECh38) positions for TF motifs in JASPAR CORE Vertebrate 2020, HOMOCOMO and Swiss Regulon. To speed up the scaning process, we divided the JASPAR2020 meme file into 12 sub-files (i.e., j1-j12.txt).
TRN/s1/s1_motif_genome_scan.sh

# Step2: get the DNase footprints of each tissue from the GTRD (the Gene Transcription Regulation Database).

## s2.1: get peak ids for each tissue.
TRN/s2/s2.1_dnase_seq_file_cellline.sh

## s2.2: combine footprint bed files for each tissue.
TRN/s2/s2.2_dnase_seq_file_cellline.sh

# Step3: Get TF-gene pairs for each tissue.

## s3.1 overlap the predicted TF binding regions (outputs of step1) with the DNase footprints of each tissue (outputs of step2).
## &
## s3.2 overlap predicted motif region of a TF with the +/-10kb region surrounding a gene's transcriptional starting site.

### use adult_cortex as an examples; for other tissues, the scripts are similar
TRN/s3/adult_cortex/tf1-12.R: correspond to predicted TF binding regions using 12 JASPAR motif subfiles.
TRN/s3/adult_cortex/tf_h1.R and tf_s1.R correspond to HOMOCOMO and Swiss Regulon motif file, seperately.

# Step4: filter TF-gene pair for each tissue.

## s4.1 get expressed genes (TPM >0.1 in at least 25% samples) in each tissue using GTEx data.
TRN/s4/s4.1_gtex_exp.R 

## s4.2 filter TF-gene pair based on expressed genes and Pearson correlation and Lasso regression.

TRN/s4/s4.2_*tissue*_tf_combine.R; 

####################################################################
# II. Construction of GCN (Gene Co-experession Networks)

WGCNA/wgcna_expgene.R

####################################################################
# III. TWAS

TWAS/TWAS_MetaXcan.sh

####################################################################
# IV. MAGMA, H-MAGMA, nMAGMA

MAGMAs/1_create_gene_and_re_annot_file.sh
MAGMAs/2_create_H-MAGMA_annot.R
MAGMAs/3_create_nMAGMA_annot.R
MAGMAs/4_running_MAGMA_H-MAGMA_nMAGMA.sh





