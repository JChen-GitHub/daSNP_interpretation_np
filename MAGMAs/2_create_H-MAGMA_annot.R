#=================================================================
#Convert all gene IDs to gene names and select all protein-coding genes
#=================================================================
setwd("your-current-working-directory")
snp_to_gene1 <- readLines("Fetal_brain.genes.annot")  # can be downloaded from https://github.com/thewonlab/H-MAGMA
library(dplyr)
gene <- c()
bp <- c()
snp <- c()
for(i in 1:length(snp_to_gene1)){
  temp <- unlist(strsplit(snp_to_gene1[i], "\t"))
  gene[i] <- temp[1]
  bp[i] <- temp[2]
  snp[i] <- paste(temp[-2:-1],collapse="\t")
}
annot0 <- data.frame(gene=gene,snp=snp,stringsAsFactors = FALSE)
#annot0=annot0[!duplicated(annot0$gene),]

library(data.table)
f1=fread("Homo_sapiens.GRCh37.87.gtf",header = F, check.names = F, sep = "\t")
f1=as.data.frame(f1)
f2=subset(f1,f1$V3=="gene")
f2=f2[,-6:-8]
chr=c("1","2","3","4","5","6","7","X","8","9","10","11","12","13","14","15","16","17","18","20","19","Y","22","21")
f2=subset(f2,f2$V1%in%chr)
library(tidyverse)
f22=separate(data = f2, col = V9, into = c("gene_id", "gene_version","gene_name","gene_source","gene_biotype"),sep = ";")
f3=f22[,c(-2,-3,-7,-9)]
f3$gene_id=gsub("gene_id \"","",f3$gene_id)
f3$gene_id=gsub("\"","",f3$gene_id) 
f3$gene_name=gsub(" gene_name \"","",f3$gene_name)
f3$gene_name=gsub("\"","",f3$gene_name)
f3$gene_biotype=gsub(" gene_biotype \"","",f3$gene_biotype)
f3$gene_biotype=gsub("\"","",f3$gene_biotype)
names(f3)=c("chr","begin","end","gene_id","gene_name","gene_biotype")

annot1=merge(annot0,f3[,c(4,5)],by.x="gene",by.y="gene_id",all.x=T)
gene_bp=read.table("protein_coding_gene.symbol.loc", header = F, check.names = F, sep = "\t")
annot1=subset(annot1,annot1$gene_name%in%gene_bp$V1)
annot1=annot1[,c(3,2)]

library(tidyr)
names(gene_bp) <- c("gene", "chr", "begin", "end")
gene_bp=tidyr::unite(gene_bp, "bp", chr, begin, sep = ":")
gene_bp=tidyr::unite(gene_bp, "bp", bp, end, sep = ":")
gene_bp=gene_bp[,c(1,2)]
annot1 <- merge(gene_bp,annot1,by.x="gene",by.y="gene_name")
annot2=annot1[!duplicated(annot1$gene),]
write.table(annot2,"Intermediate.genes.annot",row.names = F,col.names=F,sep="\t",quote=F)

#Remove duplicate snp on the same gene
data <- readLines("Intermediate.genes.annot")
data <- as.list(data)
temp <- list()
for(i in 1:length(data)){
  temp[[i]] <- unlist(strsplit(data[[i]], "\t"))
  temp[[i]] <- as.vector(temp[[i]])
  temp[[i]] <- unique(temp[[i]])
}
file.remove('Intermediate.genes.annot')
for(i in 1:length(temp)){
  cat(temp[[i]], file = "Fetal_brain_coding.genes.annot", sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
  cat("\n", file = "Fetal_brain_coding.genes.annot", sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
}
