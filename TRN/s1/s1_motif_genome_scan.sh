#!/bin/bash

# Genome scan of TF motifs

## for motifs in JASPAR 2020
fimo --text  --skip-matched-sequence ./your-working-directory/data/JASPAR/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt  your-working-directory/hg38/genome.fa > your-working-directory/JASPAR.tsv

### Here, to speed up the scaning process, we divided the JASPAR2020 meme file into 12 sub-files.
fimo --text  --skip-matched-sequence your-working-directory/data/JASPAR/j*.txt  your-working-directory/hg38/genome.fa > your-working-directory/j*.tsv


## for motifs in HOCOMOCO
fimo --text --skip-matched-sequence your-working-directory/data/GTRD/ho_addi_motif.txt  your-working-directory/hg38/genome.fa > your-working-directory/h1.tsv

## for motifs in SwissRegulon
fimo --text --skip-matched-sequence your-working-directory/data/GTRD/swiss_addi_motif.txt   your-working-directory/hg38/genome.fa > your-working-directory/s1.tsv

 
