#==============================
# Create Hi-C annotate files
#==============================
# For creating Hi-C annotate files, we assign the SNPs to genes via SNP-to-RE annotations and RE-gene regulatory pairs.
rm(list=ls())
setwd("your-current-working-directory")
snp_to_re <- readLines("1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.snp_to_re.genes.annot")
re_to_gene <- read.table("3DIV_DLPFC_RE_to_genes.txt",header = F,  sep = "\t",stringsAsFactors = F)
library(dplyr)
re <- c()
snp <- c()
for(i in 1:length(snp_to_re)){
  temp <- unlist(strsplit(snp_to_re[i], "\t"))
  re[i] <- temp[1]
  snp[i] <- paste(temp[-2:-1],collapse=" ")
}
snp_to_re2 <- data.frame(re=re,snp=snp,stringsAsFactors = FALSE)
#snp_to_re2 <- snp_to_re2[-1:-2,] 
b <- intersect(re_to_gene[,1],snp_to_re2[,1])
b <- as.vector(b)
re_to_gene <- subset(re_to_gene,re_to_gene$V1%in%b)
snp_to_re2 <- subset(snp_to_re2,snp_to_re2[,1]%in%b)
snp_to_gene <- merge(re_to_gene,snp_to_re2,by.x="V1",by.y="re",all.x=T)
snp_to_gene <- snp_to_gene[,-1]

names(snp_to_gene)=c("gene","snp")
snp_to_gene <- snp_to_gene %>% group_by(gene) %>% summarise(snp = paste(snp, collapse = " "))

# add gene coordinates
gene_bp <- read.table("protein_coding_gene.symbol.loc", header = F,  sep = "\t",stringsAsFactors = F)
names(gene_bp) <- c("gene", "chr", "begin", "end")
library(tidyr)
gene_bp <- tidyr::unite(gene_bp, "bp", chr, begin, sep = ":")
gene_bp <- tidyr::unite(gene_bp, "bp", bp, end, sep = ":")
gene_bp=gene_bp[,c(1,2)]

final_snp_to_gene <- merge(snp_to_gene,gene_bp,by.x="gene",by.y="gene")
final_snp_to_gene <- final_snp_to_gene[,c(1,3,2)]
write.table(final_snp_to_gene,"Intermediate.genes.annot",row.names = F,col.names=F,sep=" ",quote=F)

#Remove duplicate SNPs on the same gene
data <- readLines("Intermediate.genes.annot")
data <- as.list(data)
temp <- list()
for(i in 1:length(data)){
  temp[[i]] <- unlist(strsplit(data[[i]], " "))
  temp[[i]] <- as.vector(temp[[i]])
  temp[[i]] <- unique(temp[[i]])
}
file.remove('Intermediate.genes.annot')
for(i in 1:length(temp)){
  cat(temp[[i]], file = "DLPFC.Hi-C.genes.annot", sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
  cat("\n", file = "DLPFC.Hi-C.genes.annot", sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
}



#====================
#Merge four annot files
#====================
setwd("")

#Step1: Merge orig+re
#The first file
snp_to_gene1 <- readLines("1_ADHD2012_pgc.adhd.full.2012-10.hg18tohg19.ovelapg1000.snp_to_gene.genes.annot")
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
annot1 <- data.frame(gene=gene,snp=snp,stringsAsFactors = FALSE)
gene_bp1 <- data.frame(gene=gene,bp=bp,stringsAsFactors = FALSE)

#The second file
snp_to_gene2 <- readLines("DLPFC.Hi-C.genes.annot")
library(dplyr)
gene <- c()
bp <- c()
snp <- c()
for(i in 1:length(snp_to_gene2)){
  temp <- unlist(strsplit(snp_to_gene2[i], "\t"))
  gene[i] <- temp[1]
  bp[i] <- temp[2]
  snp[i] <- paste(temp[-2:-1],collapse="\t")
}
annot2 <- data.frame(gene=gene,snp=snp,stringsAsFactors = FALSE)
gene_bp2 <- data.frame(gene=gene,bp=bp,stringsAsFactors = FALSE)

#Merge
annot_comb <- rbind(annot1,annot2)
#If the annotate files contains some non-coding genes, and you want retain only coding genes
gene_bp=read.table("protein_coding_gene.symbol.loc", header = F,  sep = "\t",stringsAsFactors = F)
names(gene_bp) <- c("gene", "chr", "begin", "end")
library(tidyr)
gene_bp=tidyr::unite(gene_bp, "bp", chr, begin, sep = ":")
gene_bp=tidyr::unite(gene_bp, "bp", bp, end, sep = ":")
gene_bp=gene_bp[,c(1,2)]
#If you have retained coding genes when you generating annotation files
#gene_bp <- rbind(gene_bp1,gene_bp2)
#gene_bp <- gene_bp[!duplicated(gene_bp$gene),]

annot <- annot_comb %>% group_by(gene) %>% summarise(snp = paste(snp, collapse = "\t"))
annot <- merge(gene_bp,annot,by="gene",all=F)
write.table(annot,"Intermediate.genes.annot",row.names = F,col.names=F,sep="\t",quote=F) 

#Remove duplicate SNPs on the same gene
data <- readLines("Intermediate.genes.annot")
data <- as.list(data)
temp <- list()
for(i in 1:length(data)){
  temp[[i]] <- unlist(strsplit(data[[i]], "\t")) #\t
  temp[[i]] <- as.vector(temp[[i]])
  temp[[i]] <- unique(temp[[i]])
}
file.remove('Intermediate.genes.annot')
for(i in 1:length(temp)){
  cat(temp[[i]], file = "Cortex.orig+re.genes.annot", sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
  cat("\n", file = "Cortex.orig+re.genes.annot", sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
}


#Step2: Merge orig+re  and eqtl annot file
#The first file
snp_to_gene1 <- readLines("Cortex.orig+re.genes.annot")
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
annot1 <- data.frame(gene=gene,snp=snp,stringsAsFactors = FALSE)
gene_bp1 <- data.frame(gene=gene,bp=bp,stringsAsFactors = FALSE)

#The second file
snp_to_gene2 <- readLines("Cortex.signifpairs.finish.txt")
library(dplyr)
gene <- c()
bp <- c()
snp <- c()
for(i in 1:length(snp_to_gene2)){
  temp <- unlist(strsplit(snp_to_gene2[i], "\t"))
  gene[i] <- temp[1]
  bp[i] <- temp[2]
  snp[i] <- paste(temp[-2:-1],collapse="\t")
}
annot2 <- data.frame(gene=gene,snp=snp,stringsAsFactors = FALSE)
gene_bp2 <- data.frame(gene=gene,bp=bp,stringsAsFactors = FALSE)

#Merge
annot_comb <- rbind(annot1,annot2)
#If the annotate files contains some non-coding genes, and you want retain only coding genes
gene_bp=read.table("protein_coding_gene.symbol.loc", header = F,  sep = "\t",stringsAsFactors = F)
names(gene_bp) <- c("gene", "chr", "begin", "end")
library(tidyr)
gene_bp=tidyr::unite(gene_bp, "bp", chr, begin, sep = ":")
gene_bp=tidyr::unite(gene_bp, "bp", bp, end, sep = ":")
gene_bp=gene_bp[,c(1,2)]
#If you have retained coding genes when you generating annotation files
#gene_bp <- rbind(gene_bp1,gene_bp2)
#gene_bp <- gene_bp[!duplicated(gene_bp$gene),]

annot <- annot_comb %>% group_by(gene) %>% summarise(snp = paste(snp, collapse = "\t"))
annot <- merge(gene_bp,annot,by="gene",all=F)
write.table(annot,"Intermediate.genes.annot",row.names = F,col.names=F,sep="\t",quote=F)

#Remove duplicate SNPs on the same gene
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
  cat(temp[[i]], file = "Cortex.orig+re+eqtl.genes.annot", sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
  cat("\n", file = "Cortex.orig+re+eqtl.genes.annot", sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
}



#Step3: Merge orig+re+eqtl and wgcna annot file
#The first file
snp_to_gene1 <- readLines("Cortex.orig+re+eqtl.genes.annot")
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
annot1 <- data.frame(gene=gene,snp=snp,stringsAsFactors = FALSE)
gene_bp1 <- data.frame(gene=gene,bp=bp,stringsAsFactors = FALSE)

#The second file
snp_to_gene2 <- readLines("Cortex.wgcna.tab.genes.annot")   #can be created by codes 'Calcu_TOM_extract_signifgenepairs.R' or downloaded directly at https://github.com/sldrcyang/nMAGMA
library(dplyr)
gene <- c()
bp <- c()
snp <- c()
for(i in 1:length(snp_to_gene2)){
  temp <- unlist(strsplit(snp_to_gene2[i], "\t"))
  gene[i] <- temp[1]
  bp[i] <- temp[2]
  snp[i] <- paste(temp[-2:-1],collapse="\t")
}
annot2 <- data.frame(gene=gene,snp=snp,stringsAsFactors = FALSE)
gene_bp2 <- data.frame(gene=gene,bp=bp,stringsAsFactors = FALSE)

#Merge
annot_comb <- rbind(annot1,annot2)
#If the annotate files contains some non-coding genes, and you want retain only coding genes
gene_bp=read.table("protein_coding_gene.symbol.loc", header = F,  sep = "\t",stringsAsFactors = F)
names(gene_bp) <- c("gene", "chr", "begin", "end")
library(tidyr)
gene_bp=tidyr::unite(gene_bp, "bp", chr, begin, sep = ":")
gene_bp=tidyr::unite(gene_bp, "bp", bp, end, sep = ":")
gene_bp=gene_bp[,c(1,2)]
#If you have retained coding genes when you generating annotation files
#gene_bp <- rbind(gene_bp1,gene_bp2)
#gene_bp <- gene_bp[!duplicated(gene_bp$gene),]

annot <- annot_comb %>% group_by(gene) %>% summarise(snp = paste(snp, collapse = "\t"))
annot <- merge(gene_bp,annot,by="gene",all=F)
write.table(annot,"Intermediate.genes.annot",row.names = F,col.names=F,sep="\t",quote=F)

#Remove duplicate SNPs on the same gene
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
  cat(temp[[i]], file = "Cortex.nMAGMA.genes.annot", sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
  cat("\n", file = "Cortex.nMAGMA.genes.annot", sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
}

