library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(stringr)
library(org.Hs.eg.db)
library('glmnet')

# the predicted TF binding regions
J1 <- read.table('your-working-directory/j2.tsv',sep='\t',header = T)

# footprint file
trn_footprint <- read.table('your-working-directory/macs2_fdr/adult_brain.bed', sep='\t')
colnames(trn_footprint) <- c('chr','start','stop','id','FDR','strand')

#convert to GRanges format

j_grange <- GRanges(Rle(J1$sequence_name,lengths =rep(1,length(J1$sequence_name)) ),
                    IRanges(start=J1$start, end=J1$stop))

t_grange <- GRanges(Rle(trn_footprint$chr,lengths =rep(1,length(trn_footprint$chr)) ),
                    IRanges(start=trn_footprint$start, end=trn_footprint$stop))
#findOverlaps
tf_j <- findOverlaps(t_grange, j_grange )

# corresponding TFs and positions
tf_j <- as.data.frame(tf_j)
tf_j <- cbind(trn_footprint[tf_j$queryHits, c("chr","start","stop","strand" )], J1[tf_j$subjectHits, c('motif_alt_id') ])
colnames(tf_j)[5] <- c('motif_alt_id')


# 
tf_j_gr <- GRanges(Rle(tf_j$chr,lengths =rep(1,length(tf_j$chr)) ),
                   IRanges(start=tf_j$start, end=tf_j$stop))


# TSS +-10k
grch38 <- read.table('your-working-directory/Homo_sapiens.GRCh38.87.txt', header = T)

all_genes <- GRanges(Rle(grch38$chr,lengths =rep(1,length(grch38$chr)) ),
                   IRanges(start=grch38$start, end=grch38$end), strand=grch38$strand)

all_gene_TSS <- resize(all_genes,1)
TSS_10k <-promoters(all_gene_TSS, 10000, 10000,use.names=F)

# findOverlaps,findOverlaps with the +-10kb region surrounding a gene's transcriptional starting site.
print(4)
adult_brain_hts <- findOverlaps(tf_j_gr, TSS_10k)
adult_brain_hts <- as.data.frame(adult_brain_hts)

adult_brain_hts$queryHits <- tf_j[adult_brain_hts$queryHits,'motif_alt_id']
adult_brain_hts$subjectHits <- grch38[adult_brain_hts$subjectHits, 'symbol' ]
adult_brain_hts <- unique(adult_brain_hts)
adult_brain_hts2 <- adult_brain_hts
#brain_hts <- brain_hts_o

save(adult_brain_hts2, file='your-working-directory/adult_brain/adult_brain_hts2.RData')
