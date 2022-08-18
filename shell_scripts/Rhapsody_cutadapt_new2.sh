#!/bin/bash -i

#require cutadapt v4.1
file1=`basename $1 .fastq.gz`
file2=`basename $2 .fastq.gz`
threads=`echo $3 `
library_type=`echo $4 `
samplename=`basename $2 | sed -e 's/_R2_001.fastq.gz//g'`


#trimming adapters by cutadapt v4 from R2 reads


if [ $library_type = "WTA" ]; then

#remove ADT, R2 adapter, polyG-containing erronous R1 reads
cutadapt --cores=$threads -u -12 -U -1 -q 20 --trim-n --report=minimal -Z -n 5 \
-a "TGGCACCCGAGAATTCCA;max_error_rate=0.1" \
-a "AAGCAGTGGTATCAACGCAGA;max_error_rate=0.1" \
-g "GGGGGGGGGGGGG;min_overlap=12;max_error_rate=0.1" \
-G "XAAGCAGTGGTATCAACGCAGA;max_error_rate=0.1" -G "XG{155};min_overlap=8;max_error_rate=0.2" \
-A "A{155}X;min_overlap=12;max_error_rate=0.2" -m 43:30 \
-o ${file1}_trim.fastq.gz -p ${file2}_trim.fastq.gz \
./data/${file1}.fastq.gz ./data/${file2}.fastq.gz 1> ${samplename}_cutadapt_summary.log
sleep 1
rm ./data/${file1}.fastq.gz
rm ./data/${file2}.fastq.gz

elif [ $library_type = "targeted" ]; then

cutadapt --cores=$threads -u -12 -U -1 -q 20 --trim-n --report=minimal -Z -n 5 \
-A "A{155}X;min_overlap=12;max_error_rate=0.2" -m 43:30 \
-o ${file1}_trim.fastq.gz -p ${file2}_trim.fastq.gz \
./data/${file1}.fastq.gz ./data/${file2}.fastq.gz 1> ${samplename}_cutadapt_summary.log
sleep 1
rm ./data/${file1}.fastq.gz
rm ./data/${file2}.fastq.gz

elif [ $library_type = "hashtag" ]; then

cutadapt --cores=$threads -U -100 -q 20 --trim-n -Z --report=minimal -m 42:15 \
-o ${file1}_trim.fastq.gz -p ${file2}_trim.fastq.gz \
./data/${file1}.fastq.gz ./data/${file2}.fastq.gz 1> ${samplename}_cutadapt_summary.log

sleep 1
rm ./data/${file1}.fastq.gz
rm ./data/${file2}.fastq.gz

elif [ $library_type = "sampletag" ]; then

cutadapt --cores=$threads -u -1 -U -70 -q 20 --trim-n -Z --report=minimal -G "XGTTGTCAAGATGCTACCGTTCAGAG;min_overlap=20" -m 42:40 \
-o ${file1}_trim.fastq.gz -p ${file2}_trim.fastq.gz \
./data/${file1}.fastq.gz ./data/${file2}.fastq.gz 1> ${samplename}_cutadapt_summary.log

#cutadapt --cores=$threads -u -1 -q 20 --trim-n -Z --report=minimal -m 56:40 \
#-o ${file1}_trim.fastq.gz -p ${file2}_trim.fastq.gz \
#./data/${file1}.fastq.gz ./data/${file2}.fastq.gz 1> ${samplename}_cutadapt_summary.log

rm ./data/${file1}.fastq.gz
rm ./data/${file2}.fastq.gz


elif [ $library_type = "ADT" ]; then

##check and remove ADT adapter
cutadapt --cores=$threads -u -82 -q 20 --trim-n -Z --report=minimal \
-g ADT="XTGGCACCCGAGAATTCCA;min_overlap=17;max_error_rate=0.2" -m 15:42 \
-o {name}_${file2}.fastq.gz -p {name}_${file1}.fastq.gz  \
./data/${file2}.fastq.gz ./data/${file1}.fastq.gz 1> ${samplename}_cutadapt_summary.log

sleep 1
rm ./data/${file1}.fastq.gz
rm ./data/${file2}.fastq.gz
mv ADT_${file2}.fastq.gz ${file2}_trim.fastq.gz
mv ADT_${file1}.fastq.gz ${file1}_trim.fastq.gz
rm unknown*${file2}.fastq.gz
rm unknown*${file1}.fastq.gz

elif [ $library_type = "Abseq" ]; then

cutadapt --cores=$threads -U 12 -U -58 -q 20 --trim-n -Z --report=minimal -m 42:40 \
-o ${file1}_trim.fastq.gz -p ${file2}_trim.fastq.gz \
./data/${file1}.fastq.gz ./data/${file2}.fastq.gz 1> ${samplename}_cutadapt_summary.log

sleep 1
rm ./data/${file1}.fastq.gz
rm ./data/${file2}.fastq.gz


else
 echo "require valid library type, acceptable is WTA, targeted, hashtag or sampletag or ADT or Abseq"
 exit 1
fi



exit 0
