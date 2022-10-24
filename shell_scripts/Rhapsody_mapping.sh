#!/bin/bash -i

#require GNU parallel module, pigz, cutadapt, bowtie2, STAR v2.7.10a, Seqkit, FastQC, python3, R required packages.
#BD Rhapsody WTA proprocessing workflow for Enhanced beads with STARsolo
#Written by Shigeyuki Shichino at 20220725

###set mapping pipeline log output


#command-line option analysis by getopt (a Japanese article)
#https://atmarkit.itmedia.co.jp/ait/articles/2003/20/news002.html


echo "Check command-line parameters..."

OPTIONS=`getopt -n $(basename $0) -o 'h' -l threads:,BAM:,downsample:,dataID:,species:,library_type:,index:,in_dir:,out_dir:,expect_cells:,UMI,help -- "$@"`

#echo $OPTIONS

##set options
eval set -- "$OPTIONS"

#echo $@

# process getopt results as arguments of shell script commands.

flag_threads=false 
flag_BAM=false
flag_downsample=false
flag_dataID=false
flag_species=false
flag_library_type=false
flag_index=false
flag_in_dir=false
flag_out_dir=false
flag_help=false
flag_UMI=false
flag_expect_cells=false

while [ $# -gt 0 ]
do
   case $1 in
     --threads) flag_threads=true;
                 if [[ "$2" =~ ^[0-9]+$ ]]; then # is arguments are numeric?
                   arg_threads=$2; shift;            # to the next argument
                 else
                   flag_threads=false;  shift;       # to the next argument
                 fi;;
     --BAM) flag_BAM=true;
                 if [[ "$2" =~ ^[0-9]+$ ]]; then # is arguments are numeric?
                   arg_BAM=$2; shift;            # to the next argument
                 else 
                   flag_BAM=false; shift;            # to the next argument
                 fi;;
     --downsample) flag_downsample=true;
                 if [[ "$2" =~ ^[0-9]+$ ]]; then # is arguments are numeric?
                   arg_downsample=$2; shift;            # to the next argument
                 else 
                   flag_downsample=false;  shift;            # to the next argument
                 fi;;
     --expect_cells) flag_expect_cells=true;
                 if [[ "$2" =~ ^[0-9]+$ ]]; then # is arguments are numeric?
                   arg_expect_cells=$2; shift;            # to the next argument
                 else 
                   flag_expect_cells=false;  shift;            # to the next argument
                 fi;;
     --dataID) flag_dataID=true; dataID=$2; shift;;
     --species) flag_species=true; species=$2; shift;;
     --library_type) flag_library_type=true; library_type=$2; shift;;
     --index) flag_index=true; index=$2; shift;;
     --in_dir) flag_in_dir=true; in_dir=$2; shift;;
     --out_dir) flag_out_dir=true; out_dir=$2; shift;;
     --UMI) flag_UMI=true; shift;;
     -h | --help) flag_help=true; shift;;
     --) shift; break;;
   esac
   shift
done

#show help if specify help option or uncorrect option speficiations.
if "${flag_help}"; then
 echo "Usage: $(basename "$0") [OPTIONS...]"
 echo " -h |--help              Display help"
 echo " --threads VALUE         Number of CPU threads to use. Default is to use all threads."
 echo " --dataID VALUE          Data ID of read1 and read2 fastq.gz files. read1 (dataID_S[XX]_L00X_R1_001.fastq.gz) is cell barcode read, and read2 (dataID_S[XX]_L00X_R2_001.fastq.gz) is cDNA or tag read."
 echo " --species VALUE         Specify sample species name. Acceptable values are mmu or hsa or rat or macaca"
 echo " --library_type VALUE    Specify library type. Acceptable values are WTA, targeted, hashtag, sampletag, ADT, or Abseq."
 echo "                         set hashtag if you used BioLegend Hashtags, set sampletag if you used BD Sampletags, set ADT if you use BioLegend Totalseq antibodies or streptavidin, set Abseq if you used BD Abseq antibodies."
 echo " --expect_cells VALUE    Specify expected cell number of the library. Default=20000. Only uses WTA or targeted mode."
 echo " --index VALUE           Specify index file or folder path. Specify STAR index (for WTA data) or bowtie2 index (hashtag, sampletag, ADT, or Abseq data)." 
 echo " --BAM VALUE             BAM flag (0=do not save BAM file, 1=save BAM file, default=0)"
 echo " --in_dir VALUE          PATH to input .fastq.gz data-stored directory"
 echo " --out_dir VALUE         PATH to output directory for saving concatenated raw .fastq.gz data."
 echo " --downsample VALUE      Specify downsampling raw read number if you want to perform downsampling analysis. Default=0 (do not perform downsampling analysis)"
 echo " --UMI                   Add this option if you want to perform UMI compression by BD Rhapsody 8-base UMI. Default=OFF (no UMI compression)"
 exit 1
fi

#fuzzy correction of input/output directory name. must be terminated by "/" within pipeline.
in_dir=`echo ${in_dir}"/"`
in_dir=`echo ${in_dir} | sed -e 's/\/\//\//g'` #if backslash duplicated, correct as single backslash
out_dir=`echo ${out_dir}"/"`
out_dir=`echo ${out_dir} | sed -e 's/\/\//\//g'` #if backslash duplicated, correct as single backslash

#check input fastq file number
Nlane=`ls ${in_dir}${dataID}* | cat | wc -l`
Nlane=$((Nlane /2))

file1=`echo ${dataID}"_R1_001.fastq.gz"`
file2=`echo ${dataID}"_R2_001.fastq.gz"`
samplename=`echo ${dataID}`


hoge=./result/${samplename}_results/
if [ -d ${hoge} ]; then
sleep 1
else
mkdir ./result/${samplename}_results/
fi

LOG_OUT="./result/${samplename}_results/${samplename}_stdout.out"
exec 1> >(tee -a $LOG_OUT)

if "${flag_threads}"; then
threads=`echo $arg_threads`
threads=$((threads /2))

else
threads=`fgrep 'processor' /proc/cpuinfo | wc -l`
threads=$((threads /2))
fi

if "${flag_BAM}"; then
bam=`echo $arg_BAM`
else
bam=`echo 0`
fi

if "${flag_downsample}"; then
downsample=`echo $arg_downsample`
else
downsample=`echo 0`
fi

if "${flag_expect_cells}"; then
expect_cells=`echo $arg_expect_cells`
else
expect_cells=`echo 20000`
fi

in_dir=`echo $in_dir`
out_intermediate=`echo $out_dir`


if [ $Nlane = 0 ]; then
echo "No fastq.gz files matched with ${samplename} exist in ${in_dir} folder."
echo "Must check dataID and desired fastq.gz file names."

 echo "Usage: $(basename "$0") [OPTIONS...]"
 echo " -h |--help              Display help"
 echo " --threads VALUE         Number of CPU threads to use. Default is to use all threads."
 echo " --dataID VALUE          Data ID of read1 and read2 fastq.gz files. read1 (dataID_S[XX]_L00X_R1_001.fastq.gz) is cell barcode read, and read2 (dataID_S[XX]_L00X_R2_001.fastq.gz) is cDNA or tag read."
 echo " --species VALUE         Specify sample species name. Acceptable values are mmu or hsa or rat or macaca"
 echo " --library_type VALUE    Specify library type. Acceptable values are WTA, targeted, hashtag, sampletag, ADT, or Abseq."
 echo "                         set hashtag if you used BioLegend Hashtags, set sampletag if you used BD Sampletags, set ADT if you use BioLegend Totalseq antibodies or streptavidin, set Abseq if you used BD Abseq antibodies."
 echo " --expect_cells VALUE    Specify expected cell number of the library. Default=20000. Only uses WTA or targeted mode."
 echo " --index VALUE           Specify index file or folder path. Specify STAR index (for WTA data) or bowtie2 index (hashtag, sampletag, ADT, or Abseq data)." 
 echo " --BAM VALUE             BAM flag (0=do not save BAM file, 1=save BAM file, default=0)"
 echo " --in_dir VALUE          PATH to input .fastq.gz data-stored directory"
 echo " --out_dir VALUE         PATH to output directory for saving concatenated raw .fastq.gz data."
 echo " --downsample VALUE      Specify downsampling raw read number if you want to perform downsampling analysis. Default=0 (do not perform downsampling analysis)"
 echo " --UMI                   Add this option if you want to perform UMI compression by BD Rhapsody 8-base UMI. Default=OFF (no UMI compression)"
 
exit 1
fi

if [ $species = hsa -o $species = mmu -o $species = rat -o $species = macaca ]; then
sleep 1
else
echo "Acceptable --species options are hsa or mmu or rat or macaca."
echo "hsa(human), mmu(mouse), rat, macaca(crab-eating macaque)."
echo "Please specify correct option"


 echo "Usage: $(basename "$0") [OPTIONS...]"
 echo " -h |--help              Display help"
 echo " --threads VALUE         Number of CPU threads to use. Default is to use all threads."
 echo " --dataID VALUE          Data ID of read1 and read2 fastq.gz files. read1 (dataID_S[XX]_L00X_R1_001.fastq.gz) is cell barcode read, and read2 (dataID_S[XX]_L00X_R2_001.fastq.gz) is cDNA or tag read."
 echo " --species VALUE         Specify sample species name. Acceptable values are mmu or hsa or rat or macaca"
 echo " --library_type VALUE    Specify library type. Acceptable values are WTA, targeted, hashtag, sampletag, ADT, or Abseq."
 echo "                         set hashtag if you used BioLegend Hashtags, set sampletag if you used BD Sampletags, set ADT if you use BioLegend Totalseq antibodies or streptavidin, set Abseq if you used BD Abseq antibodies."
 echo " --expect_cells VALUE    Specify expected cell number of the library. Default=20000. Only uses WTA or targeted mode."
 echo " --index VALUE           Specify index file or folder path. Specify STAR index (for WTA data) or bowtie2 index (hashtag, sampletag, ADT, or Abseq data)." 
 echo " --BAM VALUE             BAM flag (0=do not save BAM file, 1=save BAM file, default=0)"
 echo " --in_dir VALUE          PATH to input .fastq.gz data-stored directory"
 echo " --out_dir VALUE         PATH to output directory for saving concatenated raw .fastq.gz data."
 echo " --downsample VALUE      Specify downsampling raw read number if you want to perform downsampling analysis. Default=0 (do not perform downsampling analysis)"
 echo " --UMI                   Add this option if you want to perform UMI compression by BD Rhapsody 8-base UMI. Default=OFF (no UMI compression)"
 
exit 1
fi

if [ $bam -eq 0 -o $bam -eq 1 ]; then
sleep 1
else
echo "Acceptable --BAM options are 0(no BAM output) or 1(create BAM output)."
echo "Please specify correct option"


 echo "Usage: $(basename "$0") [OPTIONS...]"
 echo " -h |--help              Display help"
 echo " --threads VALUE         Number of CPU threads to use. Default is to use all threads."
 echo " --dataID VALUE          Data ID of read1 and read2 fastq.gz files. read1 (dataID_S[XX]_L00X_R1_001.fastq.gz) is cell barcode read, and read2 (dataID_S[XX]_L00X_R2_001.fastq.gz) is cDNA or tag read."
 echo " --species VALUE         Specify sample species name. Acceptable values are mmu or hsa or rat or macaca"
 echo " --library_type VALUE    Specify library type. Acceptable values are WTA, targeted, hashtag, sampletag, ADT, or Abseq."
 echo "                         set hashtag if you used BioLegend Hashtags, set sampletag if you used BD Sampletags, set ADT if you use BioLegend Totalseq antibodies or streptavidin, set Abseq if you used BD Abseq antibodies."
 echo " --expect_cells VALUE    Specify expected cell number of the library. Default=20000. Only uses WTA or targeted mode."
 echo " --index VALUE           Specify index file or folder path. Specify STAR index (for WTA data) or bowtie2 index (hashtag, sampletag, ADT, or Abseq data)." 
 echo " --BAM VALUE             BAM flag (0=do not save BAM file, 1=save BAM file, default=0)"
 echo " --in_dir VALUE          PATH to input .fastq.gz data-stored directory"
 echo " --out_dir VALUE         PATH to output directory for saving concatenated raw .fastq.gz data."
 echo " --downsample VALUE      Specify downsampling raw read number if you want to perform downsampling analysis. Default=0 (do not perform downsampling analysis)"
 echo " --UMI                   Add this option if you want to perform UMI compression by BD Rhapsody 8-base UMI. Default=OFF (no UMI compression)"
 
exit 1
fi

if [ $library_type = WTA -o $library_type = targeted -o $library_type = hashtag -o $library_type = sampletag -o $library_type = ADT -o $library_type = Abseq ]; then
sleep 1
else
echo "Acceptable --library_type options are WTA, targeted, hashtag, sampletag, ADT, or Abseq."
echo "Please specify correct option"


 echo "Usage: $(basename "$0") [OPTIONS...]"
 echo " -h |--help              Display help"
 echo " --threads VALUE         Number of CPU threads to use. Default is to use all threads."
 echo " --dataID VALUE          Data ID of read1 and read2 fastq.gz files. read1 (dataID_S[XX]_L00X_R1_001.fastq.gz) is cell barcode read, and read2 (dataID_S[XX]_L00X_R2_001.fastq.gz) is cDNA or tag read."
 echo " --species VALUE         Specify sample species name. Acceptable values are mmu or hsa or rat or macaca"
 echo " --library_type VALUE    Specify library type. Acceptable values are WTA, targeted, hashtag, sampletag, ADT, or Abseq."
 echo "                         set hashtag if you used BioLegend Hashtags, set sampletag if you used BD Sampletags, set ADT if you use BioLegend Totalseq antibodies or streptavidin, set Abseq if you used BD Abseq antibodies."
 echo " --expect_cells VALUE    Specify expected cell number of the library. Default=20000. Only uses WTA or targeted mode."
 echo " --index VALUE           Specify index file or folder path. Specify STAR index (for WTA data) or bowtie2 index (hashtag, sampletag, ADT, or Abseq data)." 
 echo " --BAM VALUE             BAM flag (0=do not save BAM file, 1=save BAM file, default=0)"
 echo " --in_dir VALUE          PATH to input .fastq.gz data-stored directory"
 echo " --out_dir VALUE         PATH to output directory for saving concatenated raw .fastq.gz data."
 echo " --downsample VALUE      Specify downsampling raw read number if you want to perform downsampling analysis. Default=0 (do not perform downsampling analysis)"
 echo " --UMI                   Add this option if you want to perform UMI compression by BD Rhapsody 8-base UMI. Default=OFF (no UMI compression)"
 
exit 1
fi

echo "-----------------------------------------------------------------------------"


echo "command-line options for sample ${samplename} analysis:"
threads=$((threads * 2))

echo "--threads $threads"
echo "--dataID $samplename"
echo "--library_type $library_type"
echo "--BAM $bam"
echo "--downsample $downsample"
echo "--expect_cells $expect_cells"
echo "--species $species"
echo "--index $index"
echo "--in_dir $in_dir"
echo "--out_dir $out_dir"
threads=$((threads / 2))

if "${flag_UMI}"; then
echo "--UMI"
fi

echo "-----------------------------------------------------------------------------"

echo "Lane number is ${Nlane}"

echo "input read 1 files"
ls ${in_dir}${dataID}*_R1_001.fastq.gz
echo "input read 2 files"
ls ${in_dir}${dataID}*_R2_001.fastq.gz

echo "-----------------------------------------------------------------------------"



if [ $library_type = Abseq ]; then
echo "build bowtie2 index for Abseq analysis..."
bowtie2-build --threads ${threads} -q -f $index Abseq
sleep 1

mv Abseq.1.bt2 ./index/Abseq.1.bt2
mv Abseq.2.bt2 ./index/Abseq.2.bt2
mv Abseq.3.bt2 ./index/Abseq.3.bt2
mv Abseq.4.bt2 ./index/Abseq.4.bt2
mv Abseq.rev.1.bt2 ./index/Abseq.rev.1.bt2
mv Abseq.rev.2.bt2 ./index/Abseq.rev.2.bt2

index=`echo './index/Abseq'`
echo "-----------------------------------------------------------------------------"

fi

#create output directory
outname=`echo $samplename | sed -e 's/Tag//g'`
outname=`echo $outname | sed -e 's/WTA//g'`
outname=`echo $outname | sed -e 's/ADT//g'`
outname=`echo $outname | sed -e 's/Target//g'`
outname=`echo $outname | sed -e 's/TCR//g'`
outname=`echo $outname | sed -e 's/BCR//g'`
out_dir1=`echo ${out_intermediate}${outname}/`
out_dir_fastq=`echo ${out_intermediate}${outname}/fastq/`

if [ -d ${out_dir1} ]; then
sleep 1
else
mkdir ${out_dir1}
fi

if [ -d ${out_dir_fastq} ]; then
sleep 1
else
mkdir ${out_dir_fastq}
fi


start_time=`date +%s`
 echo "Start preprocessing BD Rhapsody TAS-Seq2 data..."
 echo "Copy and concat fastq files..."
 echo `date '+%y/%m/%d %H:%M:%S'`

#copy and concatenate fastq data

files="${in_dir}${dataID}*001.fastq.gz"

for filepath in ${files}

do

tmp=`basename ${filepath} .fastq.gz`

cp ${filepath} ./data/${tmp}.fastq.gz

done

ls ./data/${samplename}*_R1_001.fastq.gz | cat | xargs cat > ./data/${file1}.fastq.gz
ls ./data/${samplename}*_R2_001.fastq.gz | cat | xargs cat > ./data/${file2}.fastq.gz

mv ./data/${file1}.fastq.gz ${out_dir_fastq}${file1}.fastq.gz

mv ./data/${file2}.fastq.gz ${out_dir_fastq}${file2}.fastq.gz

echo "Finish concatenating fastq files"
echo `date '+%y/%m/%d %H:%M:%S'`

#perform downsampling fastq data if specify downsampling option.


if [ $downsample -eq 0 ]; then
sleep 1
else
echo "-----------------------------------------------------------------------------"

echo "Downsampling total reads to ${downsample} reads by seqkit..."
echo `date '+%y/%m/%d %H:%M:%S'`

tmp=`echo $((threads/Nlane*2))`
downsample1=`echo $((downsample/Nlane))`

ls ./data/${samplename}*_001.fastq.gz | cat | sort | parallel -P ${Nlane} -N2 -a \
- 'sh ./shell_scripts/Rhapsody_seqkit_downsample.sh' {1} {2} $tmp $downsample1

#remove original data

rm ./data/${samplename}*_001.fastq.gz

#change filename
files='*.part*.fastq.gz'

for filepath in $files
do
tempname1=`basename ${filepath} .fastq.gz |sed -e 's/.part_001//g'`
tempname2=`basename ${filepath} .fastq.gz |sed -e 's/_001.part_.*/_/g'`
mv $filepath ./data/${tempname1}.fastq.gz

done
echo "Finish read downsampling"
echo `date '+%y/%m/%d %H:%M:%S'`
fi



echo "-----------------------------------------------------------------------------"
#WTA or targeted annotation pipeline
if [ $library_type != "hashtag" ] && [ $library_type != "sampletag" ] && [ $library_type != "ADT" ] &&  [ $library_type != "Abseq" ]; then

#Cutadapt. correct phase shift, adapter trimming, quality filtering.
echo "Remove 5prime adapters and polyA from R2 reads..."
echo `date '+%y/%m/%d %H:%M:%S'`

tmp=`echo $((threads/Nlane*2))`
ls ./data/${samplename}*_001.fastq.gz | cat | sort | parallel -P ${Nlane} -N2 -a \
- 'sh ./shell_scripts/Rhapsody_cutadapt_new2.sh' {1} {2} $tmp $library_type
sleep 1

Rscript ./Rscripts/cutadaptlog_concatenate.R ${samplename}

end_time=`date +%s`
time_cutadapt=$((end_time - start_time))
echo "Finish read trimming"
echo `date '+%y/%m/%d %H:%M:%S'`
#splitfastq
start_time=`date +%s`
#mkdir ./temp/

echo "-----------------------------------------------------------------------------"
echo "Start base-composition analysis of sub-sequences by Seqkit and FastQC..."
echo `date '+%y/%m/%d %H:%M:%S'`

arg_read1=`ls ${dataID}*_R1_001_trim.fastq.gz | cat | sort | sed -n 1p `
arg_read2=`ls ${dataID}*_R2_001_trim.fastq.gz | cat | sort | sed -n 1p `

seqkit sample -p 0.1 $arg_read1 | seqkit head -n 5000000 > ${samplename}_fastqc_R1.fastq
seqkit sample -p 0.1 $arg_read2 | seqkit head -n 5000000 > ${samplename}_fastqc_R2.fastq
pigz -p 16 ${samplename}_fastqc_R1.fastq
pigz -p 16 ${samplename}_fastqc_R2.fastq

sleep 1


#check fastq quality by FastQC by using subset of fastq data
fastqc -t 2 --nogroup -q -o ./data ${samplename}_fastqc_R1.fastq.gz ${samplename}_fastqc_R2.fastq.gz
sleep 1
unzip -qq ./data/${samplename}_fastqc_R1_fastqc.zip
unzip -qq ./data/${samplename}_fastqc_R2_fastqc.zip
mv ${samplename}_fastqc_R1_fastqc/Images/per_base_sequence_content.png ./data/${samplename}_per_base_sequence_content_R1.png
mv ${samplename}_fastqc_R2_fastqc/Images/per_base_sequence_content.png ./data/${samplename}_per_base_sequence_content_R2.png
rm ./data/${samplename}_fastqc_R1_fastqc.zip
rm ./data/${samplename}_fastqc_R2_fastqc.zip
rm ./data/${samplename}_fastqc_R1_fastqc.html
rm ./data/${samplename}_fastqc_R2_fastqc.html
rm -Rf ${samplename}_fastqc_R1_fastqc
rm -Rf ${samplename}_fastqc_R2_fastqc

#remove chunked data
rm ${samplename}_fastqc_R1.fastq.gz
rm ${samplename}_fastqc_R2.fastq.gz


end_time=`date +%s`
time_splitfastq=$((end_time - start_time))
echo "Finish FastQC analysis"
echo `date '+%y/%m/%d %H:%M:%S'`
#mapping by STARsolo
echo "-----------------------------------------------------------------------------"
start_time=`date +%s`
echo "Start STARsolo analysis..."
echo `date '+%y/%m/%d %H:%M:%S'`
#make directory
hoge=`echo "./result/${samplename}_results"`
if [ -d ${hoge} ]; then
rm -fR $hoge
fi

mkdir ./result/${samplename}_results
mkdir ./result/${samplename}_results/matrix
mkdir ./result/${samplename}_results/processed
mkdir ./result/${samplename}_results/stats
mkdir ./result/${samplename}_results/plots
mkdir ./result/${samplename}_results/Gene
mkdir ./result/${samplename}_results/velocyto

tmp1=`echo $((threads/Nlane*2))`


ls ${samplename}*_001_trim.fastq.gz | cat | sort | parallel -P ${Nlane} -N2 -a \
- 'sh ./shell_scripts/Rhapsody_STARsolo.sh' {1} {2} ${tmp1} ${index} ${flag_UMI} ${bam} ${samplename}


echo "End STARsolo analysis..."
end_time=`date +%s`
time_annotate=$((end_time - start_time))

echo "-----------------------------------------------------------------------------"

echo "generate gene-cells martix and identify valid cell barcodes..."
start_time=`date +%s`
sleep 1

#concatenate STARsolo count table and Summary statistics

ls -d ${samplename}*L00* | cat > dir_list.txt

Rscript ./Rscripts/STARsolo_concatenate.R ${samplename} ${Nlane}

rm -fR dir_list.txt

if [ $bam = 1 ]; then

echo "concatenate bam files..."
ls ${samplename}*.bam | cat | xargs samtools merge -@${threads} ${samplename}.BAM
mv ${samplename}.BAM ./result/${samplename}_results/processed/${samplename}.bam
rm ${samplename}*.bam


sleep 1
fi

#create html report by Rmarkdown package


if "${flag_UMI}"; then
flag_UMI=TRUE
else
flag_UMI=FALSE
fi

if [ $library_type = WTA ]; then

echo "Perform dropkick analysis..."
#perform dropkick analysis to estimate background beads and valid cells
#https://github.com/KenLauLab/dropkick
#https://genome.cshlp.org/content/31/10/1742

dropkick run -j 32 --seed 656565 --csv --n-ambient 20 --min-genes 200 ./result/${samplename}_results/processed/${samplename}_counts_top80000BC.h5ad
mv ${samplename}_counts_top80000BC_score.png ./result/${samplename}_results/plots/${samplename}_counts_top80000BC_score.png
mv ${samplename}_counts_top80000BC_dropkick_scores.csv ./result/${samplename}_results/processed/${samplename}_counts_top80000BC_dropkick_scores.csv
pigz -p 16 ./result/${samplename}_results/processed/${samplename}_counts_top80000BC_dropkick_scores.csv

dropkick qc --n-ambient 20 ./result/${samplename}_results/processed/${samplename}_counts_top80000BC.h5ad
mv ${samplename}_counts_top80000BC_qc.png ./result/${samplename}_results/plots/${samplename}_counts_top80000BC_qc.png
rm ${samplename}_counts_top80000BC_coef.png
rm ${samplename}_counts_top80000BC_dropkick.h5ad

pigz -p 16 ./result/${samplename}_results/processed/${samplename}_counts_top80000BC.h5ad

Rscript -e "rmarkdown::render('./Rscripts/BDWTA_matrixCreate_STARsolo.Rmd', output_file = '${samplename}_mapping_report.html', output_dir='./result/${samplename}_results/', params=list(samplename='${samplename}', species='${species}', umi='${flag_UMI}', expectCells='${expect_cells}'))"
rm -R ./result/${samplename}_results/${samplename}_mapping_report_files
rm -R ./result/${samplename}_results/Gene/

else

Rscript -e "rmarkdown::render('./Rscripts/BDWTA_matrixCreate_STARsolo_targeted.Rmd', output_file = '${samplename}_mapping_report.html', output_dir='./result/${samplename}_results/', params=list(samplename='${samplename}', species='${species}', umi='${flag_UMI}', expectCells='${expect_cells}'))"
rm -R ./result/${samplename}_results/${samplename}_mapping_report_files
rm -R ./result/${samplename}_results/Gene/

fi

sleep 1
rm ./result/*${samplename}*.log

#move processed files
mv ./data/${samplename}_per_base_sequence_content_R1.png ./result/${samplename}_results/plots/${samplename}_per_base_sequence_content_R1.png
mv ./data/${samplename}_per_base_sequence_content_R2.png ./result/${samplename}_results/plots/${samplename}_per_base_sequence_content_R2.png

end_time=`date +%s`
time_reshaping=$((end_time - start_time))


rm -fR ${samplename}*L00*

#create anndata for scVelo analysis
if [ $library_type = WTA ]; then

echo "create anndata for scVelo analysis..."
python3 ./Rhapsody_python_new/RhapsodyPython/apps/STARsoloToAdata.py --resultPath ./result/${samplename}_results/
fi


echo "Preprocessing finished." 
echo "Cutadapt trimming taken ${time_cutadapt} seconds" 
echo "FastQC analysis taken ${time_splitfastq} seconds" 
echo "STARsolo analysis taken ${time_annotate} seconds" 
echo "Generating genes-cells matrix taken ${time_reshaping} seconds"
time=$((time_cutadapt + time_splitfastq + time_annotate + time_reshaping))
echo "Total running time is ${time} seconds"
echo `date '+%y/%m/%d %H:%M:%S'`

#cp ./result/${samplename}_results/${samplename}_stdout.log ./result/${samplename}_results/stats/${samplename}_stdout.log


else


##---------------------------------------------------------------------------------------------
#tag mapping pipeline


echo "Remove 5prime Adapters from R2 reads..."

tmp=`echo $((threads/Nlane*2))`
ls ./data/${samplename}*_001.fastq.gz | cat | sort | parallel -P ${Nlane} -N2 -a \
- 'sh ./shell_scripts/Rhapsody_cutadapt_new2.sh' {1} {2} $tmp $library_type
sleep 1

Rscript ./Rscripts/cutadaptlog_concatenate.R ${samplename}

end_time=`date +%s`
time_cutadapt=$((end_time - start_time))

echo "-----------------------------------------------------------------------------"
echo "Start base-composition analysis of sub-sequences by Seqkit and FastQC..."
echo `date '+%y/%m/%d %H:%M:%S'`

arg_read1=`ls ${dataID}*_R1_001_trim.fastq.gz | cat | sort | sed -n 1p `
arg_read2=`ls ${dataID}*_R2_001_trim.fastq.gz | cat | sort | sed -n 1p `

seqkit sample -p 0.1 $arg_read1 | seqkit head -n 1000000 > ${samplename}_fastqc_R1.fastq
seqkit sample -p 0.1 $arg_read2 | seqkit head -n 1000000 > ${samplename}_fastqc_R2.fastq
pigz -p 16 ${samplename}_fastqc_R1.fastq
pigz -p 16 ${samplename}_fastqc_R2.fastq

sleep 1


#check fastq quality by FastQC by using subset of fastq data
fastqc -t 2 --nogroup -q -o ./data ${samplename}_fastqc_R1.fastq.gz ${samplename}_fastqc_R2.fastq.gz
sleep 1
unzip -qq ./data/${samplename}_fastqc_R1_fastqc.zip
unzip -qq ./data/${samplename}_fastqc_R2_fastqc.zip
mv ${samplename}_fastqc_R1_fastqc/Images/per_base_sequence_content.png ./data/${samplename}_per_base_sequence_content_R1.png
mv ${samplename}_fastqc_R2_fastqc/Images/per_base_sequence_content.png ./data/${samplename}_per_base_sequence_content_R2.png
rm ./data/${samplename}_fastqc_R1_fastqc.zip
rm ./data/${samplename}_fastqc_R2_fastqc.zip
rm ./data/${samplename}_fastqc_R1_fastqc.html
rm ./data/${samplename}_fastqc_R2_fastqc.html
rm -Rf ${samplename}_fastqc_R1_fastqc
rm -Rf ${samplename}_fastqc_R2_fastqc

#remove chunked data
rm ${samplename}_fastqc_R1.fastq.gz
rm ${samplename}_fastqc_R2.fastq.gz

######--------------------------------------------------------------------

#splitfastq
start_time=`date +%s`
#mkdir ./temp/
echo "-----------------------------------------------------------------------------"


echo "Start Splitting Fastq into ${threads} files by Seqkit..."
echo `date '+%y/%m/%d %H:%M:%S'`

tmp=`echo $((threads/Nlane*2))`
ls ${samplename}*_001_trim.fastq.gz | cat | sort | parallel -P ${Nlane} -N2 -a \
- 'sh ./shell_scripts/Rhapsody_seqkit.sh' {1} {2} $tmp
sleep 1

rm ${samplename}*_trim.fastq.gz

sleep 1

#change filename
files='*.part*.fastq.gz'

for filepath in $files
do
tempname1=`basename ${filepath} .fastq.gz |sed -e 's/.*_trim.part_//g'`
tempname2=`basename ${filepath} .fastq.gz |sed -e 's/_001_trim.part_.*/_/g'`
mv $filepath ${tempname1}_${tempname2}.fastq.gz

done

sleep 1
end_time=`date +%s`
time_splitfastq=$((end_time - start_time))


#annotateR1 in parallel
echo "-----------------------------------------------------------------------------"

start_time=`date +%s`
echo "Start R1 Annotation..."
echo `date '+%y/%m/%d %H:%M:%S'`

tmp=`echo $((threads/Nlane*2))`
ls *${samplename}*_R1_.fastq.gz | cat | sort | parallel -P ${tmp} -a \
- 'python3 ./Rhapsody_python_new/RhapsodyPython/apps/AnnotateR1.py' --R1 {} --cutadapt_log ./result/${samplename}_cutadapt.log
sleep 1

rm *_R1_.fastq.gz
rm AnnotateR1.log 

ls *${samplename}*.csv.gz | cat | sort | parallel -P ${threads} -a \
- 'Rscript ./Rscripts/R1_seq_retrieve.R' {}
sleep 1

echo "End R1 Annotation..."
end_time=`date +%s`
time_annotateR1=$((end_time - start_time))

#AnnotateR2 and AddtoSam
echo "-----------------------------------------------------------------------------"

start_time=`date +%s`
echo "Start R2 Annotation and add R1 information to mapped R2 files..."
echo `date '+%y/%m/%d %H:%M:%S'`

tmp=`echo $((threads/Nlane*2))`

ls *${samplename}*_R2_.fastq.gz | cat | sort | parallel -P ${tmp} -a \
- 'sh ./shell_scripts/Rhapsody_bowtie2_addtoBam_new.sh' {} $species $bam $library_type $index
sleep 1

rm *_R1_error_count_table.npy
rm *_R1_read_count_breakdown.json

Rscript ./Rscripts/mappinglog_concatenate.R ${samplename}
sleep 1
rm *_samtools.log

echo "End R1 and R2 Annotation"

echo "concatenate annotated files..."
cat *_count_R2.txt | gawk '{sum[$1"\t"$2]+=$3} END {for (name in sum) print name"\t"sum[name]}'> ./result/${samplename}_final_count.txt

sleep 1
pigz -f ./result/${samplename}_final_count.txt
rm *_R2.txt
end_time=`date +%s`
time_annotateR2_and_AddtoSam=$((end_time - start_time))
echo "-----------------------------------------------------------------------------"

echo "generate gene-cells martix tables..."
start_time=`date +%s`
sleep 1

#summarize count table (no UMI compression)
#create html report by Rmarkdown package

Rscript -e "rmarkdown::render('./Rscripts/BDTag_matrixCreate_new1.Rmd', output_file = '${samplename}_mapping_report.html', output_dir='./result/${samplename}_results/', params=list(samplename='${samplename}', species='${species}', adapter='${library_type}'))"

rm -R ./result/${samplename}_results/${samplename}_mapping_report_files
sleep 1
rm ./result/*${samplename}*.log

#move processed files
mv ./result/${samplename}_final_count.txt.gz ./result/${samplename}_results/processed/${samplename}_final_count.txt.gz
mv ./data/${samplename}_per_base_sequence_content_R1.png ./result/${samplename}_results/plots/${samplename}_per_base_sequence_content_R1.png
mv ./data/${samplename}_per_base_sequence_content_R2.png ./result/${samplename}_results/plots/${samplename}_per_base_sequence_content_R2.png

end_time=`date +%s`
time_reshaping=$((end_time - start_time))


echo "Preprocessing finished." 
echo "Cutadapt ${time_cutadapt} seconds" 
echo "splitfastq ${time_splitfastq} seconds" 
echo "annotateR1 ${time_annotateR1} seconds" 
echo "annotateR2_and_AddtoSam ${time_annotateR2_and_AddtoSam} seconds" 
echo "generating matrix ${time_reshaping} seconds" 
time=$((time_cutadapt + time_splitfastq + time_annotateR1 + time_annotateR2_and_AddtoSam + time_reshaping))
echo "Total running time is ${time} seconds"
echo `date '+%y/%m/%d %H:%M:%S'`


fi
echo "--------------------------------------------------------"


exit 0
