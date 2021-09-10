#!bin/bash

######################
while read Line
do
#echo $Line
#ls your-working-directory/macs2_fdr/${Line}*20.bed|wc -l
cat your-working-directory/macs2_fdr/${Line}*20.bed  >>your-working-directory/macs2_fdr/hippocampus.bed

done <your-working-directory/data/GTRD/hippocampus.txt

######################
while read Line
do
#echo $Line
#ls your-working-directory/macs2_fdr/${Line}*20.bed|wc -l
cat your-working-directory/macs2_fdr/${Line}*20.bed  >>your-working-directory/macs2_fdr/intestine.bed

done <your-working-directory/data/GTRD/intestine.txt


######################
while read Line
do
#echo $Line
cat your-working-directory/macs2_fdr/${Line}*20.bed  >>your-working-directory/macs2_fdr/colon.bed

done <your-working-directory/data/GTRD/colon.txt

# ######################
while read Line
do
#echo $Line
cat your-working-directory/macs2_fdr/${Line}*20.bed  >>your-working-directory/macs2_fdr/lung.bed

done <your-working-directory/data/GTRD/lung.txt

#######################
while read Line
do
#echo $Line
cat your-working-directory/macs2_fdr/${Line}*20.bed  >>your-working-directory/macs2_fdr/liver.bed

done <your-working-directory/data/GTRD/liver.txt


#####################################

while read Line
do
#echo $Line
cat your-working-directory/macs2_fdr/${Line}*20.bed  >>your-working-directory/macs2_fdr/stomach.bed

done <your-working-directory/data/GTRD/stomach.txt

######################
while read Line
do
#echo $Line
#ls your-working-directory/macs2_fdr/${Line}*20.bed|wc -l
cat your-working-directory/macs2_fdr/${Line}*20.bed  >>your-working-directory/macs2_fdr/brain.bed

done <your-working-directory/data/GTRD/brain.txt

