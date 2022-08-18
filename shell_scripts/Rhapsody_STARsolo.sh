#!/bin/bash -i


file1=`basename $1 _trim.fastq.gz`
file2=`basename $2 _trim.fastq.gz`
threads=`echo $3 `
index=$4
samplename=`basename $2 | sed -e 's/_R2_001_trim.fastq.gz//g' `

flag_UMI=$5
bam=$6
samplename2=$7


if "${flag_UMI}"; then

if [ $bam = 0 ]; then
STAR --runThreadN ${threads} --genomeDir ${index} \
--readFilesIn ${file2}_trim.fastq.gz ${file1}_trim.fastq.gz \
--readFilesCommand gunzip -c \
--clipAdapterType CellRanger4 \
--outFileNamePrefix ./${samplename}/ \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
--outFilterMultimapScoreRange 0 --seedSearchStartLmax 30 \
--soloCellFilter None \
--soloOutFileNames ./${samplename}/${samplename}_ \
--soloUMIdedup Exact \
--soloMultiMappers Rescue \
--soloFeatures Gene GeneFull \
--soloAdapterSequence NNNNNNNNNGTGANNNNNNNNNGACA \
--soloCBmatchWLtype EditDist_2 \
--soloCBwhitelist "./reference/BD_CB/BD_CLS1.txt" "./reference/BD_CB/BD_CLS2.txt" "./reference/BD_CB/BD_CLS3.txt" \
--soloType CB_UMI_Complex \
--soloUMIlen 8 \
--soloCBposition 2_0_2_8 2_13_2_21 3_1_3_9 \
--soloUMIposition 3_10_3_17 
else
STAR --runThreadN ${threads} --genomeDir ${index} \
--readFilesIn ${file2}_trim.fastq.gz ${file1}_trim.fastq.gz \
--readFilesCommand gunzip -c \
--clipAdapterType CellRanger4 \
--outFileNamePrefix ./${samplename}/ \
--outSAMtype BAM SortedByCoordinate \
--outBAMsortingThreadN ${threads} \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--outSAMunmapped Within \
--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
--outFilterMultimapScoreRange 0 --seedSearchStartLmax 30 \
--soloCellFilter None \
--soloOutFileNames ./${samplename}/${samplename}_ \
--soloUMIdedup Exact \
--soloMultiMappers Rescue \
--soloFeatures Gene GeneFull \
--soloAdapterSequence NNNNNNNNNGTGANNNNNNNNNGACA \
--soloCBmatchWLtype EditDist_2 \
--soloCBwhitelist "./reference/BD_CB/BD_CLS1.txt" "./reference/BD_CB/BD_CLS2.txt" "./reference/BD_CB/BD_CLS3.txt" \
--soloType CB_UMI_Complex \
--soloUMIlen 8 \
--soloCBposition 2_0_2_8 2_13_2_21 3_1_3_9 \
--soloUMIposition 3_10_3_17 
fi

else

if [ $bam = 0 ]; then
STAR --runThreadN ${threads} --genomeDir ${index} \
--readFilesIn ${file2}_trim.fastq.gz ${file1}_trim.fastq.gz \
--readFilesCommand gunzip -c \
--clipAdapterType CellRanger4 \
--outFileNamePrefix ./${samplename}/ \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
--outFilterMultimapScoreRange 0 --seedSearchStartLmax 30 \
--soloCellFilter None \
--soloOutFileNames ./${samplename}/${samplename}_ \
--soloUMIdedup NoDedup \
--soloMultiMappers Rescue \
--soloFeatures Gene GeneFull \
--soloAdapterSequence NNNNNNNNNGTGANNNNNNNNNGACA \
--soloCBmatchWLtype EditDist_2 \
--soloCBwhitelist "./reference/BD_CB/BD_CLS1.txt" "./reference/BD_CB/BD_CLS2.txt" "./reference/BD_CB/BD_CLS3.txt" \
--soloType CB_UMI_Complex \
--soloUMIlen 8 \
--soloCBposition 2_0_2_8 2_13_2_21 3_1_3_9 \
--soloUMIposition 3_10_3_17 
else
STAR --runThreadN ${threads} --genomeDir ${index} \
--readFilesIn ${file2}_trim.fastq.gz ${file1}_trim.fastq.gz \
--readFilesCommand gunzip -c \
--clipAdapterType CellRanger4 \
--outFileNamePrefix ./${samplename}/ \
--outSAMtype BAM SortedByCoordinate \
--outBAMsortingThreadN ${threads} \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--outSAMunmapped Within \
--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
--outFilterMultimapScoreRange 0 --seedSearchStartLmax 30 \
--soloCellFilter None \
--soloOutFileNames ./${samplename}/${samplename}_ \
--soloUMIdedup NoDedup \
--soloMultiMappers Rescue \
--soloFeatures Gene GeneFull \
--soloAdapterSequence NNNNNNNNNGTGANNNNNNNNNGACA \
--soloCBmatchWLtype EditDist_2 \
--soloCBwhitelist "./reference/BD_CB/BD_CLS1.txt" "./reference/BD_CB/BD_CLS2.txt" "./reference/BD_CB/BD_CLS3.txt" \
--soloType CB_UMI_Complex \
--soloUMIlen 8 \
--soloCBposition 2_0_2_8 2_13_2_21 3_1_3_9 \
--soloUMIposition 3_10_3_17 
fi

fi


rm ./${samplename}/*.out
rm ./${samplename}/*.tab
rm ${file1}_trim.fastq.gz
rm ${file2}_trim.fastq.gz

sleep 1

if [ $bam = 1 ]; then 
mv ./${samplename}/Aligned.sortedByCoord.out.bam ${samplename}.bam
else
rm ./${samplename}/*.bam
fi

pigz -p 8 -f ./${samplename}/${samplename}/${samplename}_GeneFull/raw/*
mv ./${samplename}/${samplename}/${samplename}_GeneFull/raw/UniqueAndMult-Rescue.mtx.gz ./${samplename}/${samplename}/${samplename}_GeneFull/raw/matrix.mtx.gz
mv ./${samplename}/${samplename}/${samplename}_GeneFull/Summary.csv ./result/${samplename2}_results/stats/${samplename}_Summary.csv

pigz -p 8 -f ./${samplename}/${samplename}/${samplename}_Gene/raw/*
mv ./${samplename}/${samplename}/${samplename}_Gene/raw/UniqueAndMult-Rescue.mtx.gz ./${samplename}/${samplename}/${samplename}_Gene/raw/matrix.mtx.gz

exit 0

