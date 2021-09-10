#get peak ids for each tissue.

gtrd <- read.table('your-working-directory/data/GTRD/METADATA.txt',sep='\t',header = T,quote = '',stringsAsFactors = F)
GTRD_human <- subset(gtrd,species=='Homo sapiens')


fetal_brain_cortex <- subset(GTRD_human,  cell_type=='fetal brain')[,2]
adult_brain_cortex <- subset(GTRD_human,  cell_type=='frontal cortex')[,2]
liver <- subset(GTRD_human, cell_type%in%c('liver','right lobe of liver') )[,2]
lung <- subset(GTRD_human, cell_type%in%c('lung','upper lobe of left lung'))[,2]
sto <- subset(GTRD_human, cell_type%in%c('stomach'))[,2]
intestine <- subset(GTRD_human,  cell_type=='small intestine')[,2]
hippocampus <- subset(GTRD_human,  cell_type%in%c('astrocytes of the hippocampus','HA-h (astrocytes-hippocampal)'))[,2]

write.table(hippocampus, file='your-working-directory/data/GTRD/hippocampus.txt',quote = F, row.names = F,col.names=F)
write.table(intestine, file='your-working-directory/data/GTRD/intestine.txt',quote = F, row.names = F,col.names=F)
write.table(adult_brain_cortex, file='your-working-directory/data/GTRD/adult_brain_cortex.txt',quote = F, row.names = F,col.names=F)
write.table(fetal_brain_cortex, file='your-working-directory/data/GTRD/fetal_brain_cortex.txt',quote = F, row.names = F,col.names=F)
write.table(liver, file='your-working-directory/data/GTRD/liver.txt',quote = F, row.names = F,col.names=F)
write.table(lung, file='your-working-directory/data/GTRD/lung.txt',quote = F, row.names = F,col.names=F)
write.table(sto, file='your-working-directory/data/GTRD/stomach.txt',quote = F, row.names = F,col.names=F)



