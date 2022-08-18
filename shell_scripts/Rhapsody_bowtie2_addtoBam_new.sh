#!/bin/bash -i

CMDNAME=`basename $0`
if [ $# -ne 5 ]; then
 echo "Usage: sh $CMDNAME Read2.fastq.gz <species> <BAM flag>"
 echo "read2 is mRNA reads" 
 echo "example : sh $CMDNAME day00-1_S1_R2_001.fastq.gz hsa 0 <adapter> <index>" 1>&5
 exit 1
fi

file1=`basename $1 _R2_.fastq.gz`
species=`echo $2`
adapter=`echo $4`
index=`echo $5`

#bowtie2
if [ $adapter != "hashtag" ] && [ $adapter != "sampletag" ] && [ $adapter != "ADT" ]; then
 bowtie2 -p 2 -N 1 --very-sensitive-local --norc --seed 656565 --reorder -x ${index}\
 -U $1|samtools view -@2 -Shb - > ${file1}_mapping_R2.BAM
 sleep 1

 samtools flagstat -@2 ${file1}_mapping_R2.BAM | sed '2,4d' - | sed '3,10d' | sed -e 's/ .*//g' > ${file1}_samtools.log
 sleep 1
 rm $1

elif [ $adapter = "hashtag" ]; then
bowtie2 -p 2 -D 50 -R 20 -N 0 -L 8 -i S,1,0.75 --norc --seed 656565 --reorder -x ${index} --trim-to 3:21\
 --score-min L,-9,0 --mp 3,3 --np 3 --rdg 3,3 -U $1|samtools view -@ 2 -Shb - > ${file1}_mapping_R2.BAM
 sleep 1
  
 samtools flagstat -@2 ${file1}_mapping_R2.BAM | sed '2,4d' - | sed '3,10d' | sed -e 's/ .*//g' > ${file1}_samtools.log
 sleep 1
 rm $1

elif [ $adapter = "sampletag" ]; then
bowtie2 -p 2 -D 20 -R 3 -N 0 -L 14 -i S,1,0.75 --norc --seed 656565 --reorder -x ${index} --trim-to 3:40\
 --score-min L,-9,0 --mp 3,3 --np 3 --rdg 3,3 -U $1|samtools view -@ 2 -Shb - > ${file1}_mapping_R2.BAM
 sleep 1
  
 samtools flagstat -@2 ${file1}_mapping_R2.BAM| sed '2,4d' - | sed '3,10d' | sed -e 's/ .*//g' > ${file1}_samtools.log
 sleep 1
 rm $1

elif [ $adapter = "ADT" ]; then
bowtie2 -p 2 -D 50 -R 20 -N 0 -L 8 -i S,1,0.75 --norc --seed 656565 --reorder -x ${index} --trim-to 3:21\
 --score-min L,-9,0 --mp 3,3 --np 3 --rdg 3,3 -U $1|samtools view -@ 2 -Shb - > ${file1}_mapping_R2.BAM
 sleep 1
  
 samtools flagstat -@2 ${file1}_mapping_R2.BAM | sed '2,4d' - | sed '3,10d' | sed -e 's/ .*//g' > ${file1}_samtools.log
 sleep 1
 rm $1
 
 elif [ $adapter = "Abseq" ]; then
bowtie2 -p 2 -D 20 -R 3 -N 0 -L 8 -i S,1,0.75 --norc --seed 656565 --reorder -x ./index/Abseq --trim-to 3:39\
 --score-min L,-9,0 --mp 3,3 --np 3 --rdg 3,3 -U $1|samtools view -@ 2 -Shb - > ${file1}_mapping_R2.BAM
 sleep 1
  
 samtools flagstat -@2 ${file1}_mapping_R2.BAM | sed '2,4d' - | sed '3,10d' | sed -e 's/ .*//g' > ${file1}_samtools.log
 sleep 1
 rm $1

fi

if [ $3 = 1 ]; then
#retain mapping info, cell BC and MI-containing BAM files 
python3 ./Rhapsody_python_new/RhapsodyPython/apps/AddtoBam_gzip.py --annot-R1 ${file1}_Annotation_R1.csv.gz --bam ${file1}_mapping_R2.BAM 


sleep 1

samtools view -@2 -h ${file1}_Annotated_mapping_R2.BAM | grep -v -e '^@\|^#\|^$\' - | grep -v 'CB:Z:0' | mawk '{print $3"\t"$(NF-1)"\t"$NF}' | grep -v '^\*' | sed -e 's/^\([^_]*\)/\1\t/g' -e 's/CB:Z://g' -e 's/MR:Z://g' | sed -e 's/\t_/\t/g' | sed -e 's/_[^\t]*//g' | mawk '{count[$0]++} END {for (name in count) print name"\t"count[name]}' > gathered_${file1}_count_R2.txt 


elif [ $3 = 0 ]; then
#parse BAM files directly to txt with pybam module
python3 ./Rhapsody_python_new/RhapsodyPython/apps/AddtoBam_pysam_gzip.py --annot-R1 ${file1}_Annotation_R1.csv.gz --bam ${file1}_mapping_R2.BAM 

sleep 1

unpigz -p 1 -dc ${file1}_Annotated_mapping_R2.txt.gz | grep -v '\*' - | sed -e 's/^\([^_]*\)/\1\t/g' | sed -e 's/\t_/\t/g' | sed -e 's/\(\t[^_]*_[^_]*\)_[^\t]*/\t\1/g' | mawk '{count[$1"\t"$3]++} END {for (name in count) print name"\t"count[name]}' > gathered_${file1}_count_R2.txt

sleep 1
rm ${file1}_Annotated_mapping_R2.txt.gz

else
sleep 1
fi

exit 0
