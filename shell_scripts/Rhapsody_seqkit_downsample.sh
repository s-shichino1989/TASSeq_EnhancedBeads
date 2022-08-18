#!/bin/bash -i

#require GNU parallel module and seqkit
#wrapper of splitting fastq files

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
seqkit split2 -j $3 -s $4 -f --quiet -O ${WORKING_DIR}/${temp}/ -1 $1 -2 $2 2> /dev/null

sleep 1


files=`echo ${WORKING_DIR}/${temp}/*.fastq.gz`

WORKING_DIR=$(dirname $(cd $(dirname $0); pwd))
#cho ${WORKING_DIR}  ; 

ls ${WORKING_DIR}/${temp}/*.part_001.fastq.gz | cat | xargs -I% mv % ${WORKING_DIR}

#for filepath in $files
#do

#mv $filepath ./ 
#echo $filepath ${WORKING_DIR}
#mv $filepath ${WORKING_DIR}

#done

sleep 1

rm -fR ${WORKING_DIR}/${temp}/
exit 0

#end of the script
