#!/bin/bash -i

#require GNU parallel module and seqkit
#wrapper of splitting fastq files

CMDNAME=`basename $0`
if [ $# -ne 3 ]; then
 echo "Usage: sh $CMDNAME Read1.fastq.gz Read2.fastq.gz t<num_threads> <adapter>"
 echo "fastq.gz files are stored in ./data folder, read1 is cell barcode reads and read2 is mRNA reads" 
 echo "example : sh $CMDNAME day00-1_S1_L001_R1_001.fastq.gz day00-1_S1_L001_R2_001.fastq.gz 16" 1>&3
 exit 1
fi

file1=`basename $1 .fastq.gz`
file2=`basename $2 .fastq.gz`
temp=`basename $1 .fastq.gz | sed -e 's/_R1_001//g'`

#load required modules

echo "make temporary new folder: ${temp}"

#mkdir ./${temp} &
if [ ! -d ${temp} ]; then
     mkdir ${temp}
else
    rm -rf  ${temp}
fi

WORKING_DIR=$(dirname $(cd $(dirname $0); pwd))

sleep 1

#splitting fastq files by seqkit
seqkit split2 -j $3 -p $3 -f --quiet -O ${WORKING_DIR}/${temp}/ -1 ${file1}.fastq.gz -2 ${file2}.fastq.gz 2> /dev/null

sleep 1


files=`echo ${WORKING_DIR}/${temp}/*.fastq.gz`

WORKING_DIR=$(dirname $(cd $(dirname $0); pwd))
#cho ${WORKING_DIR}  ; 

ls ${WORKING_DIR}/${temp}/*.fastq.gz | cat | xargs -I% mv % ${WORKING_DIR}

#for filepath in $files
#do

#mv $filepath ./ 
#echo $filepath ${WORKING_DIR}
#mv $filepath ${WORKING_DIR}

#done

sleep 1

rm -R ${WORKING_DIR}//${temp}/
exit 0

#end of the script
