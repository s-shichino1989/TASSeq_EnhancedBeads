# TAS-Seq mapping pipeline
Analytical package for TAS-Seq / BD Rhapsody WTA data analysis (BD Rhapsody enhanced beads)

Clone repository, extract, and change folder name as "Rhapsody_analysis".

Please first see setup.sh for creating analytical environment.
This package is tested for Ubuntu 20.04 LTS.

CAUTION!! R environment is forced to clean and re-install R 4.2.1.
You could setup analytical environment by moving on to the Rhapsody_analysis directory and typing 

sudo sh setup.sh

if you do not want to remove your current R environment, install required packages manually 
by seeing setup.sh and ./Rscripts/setup.R


Once steup.sh is finished, add permission to main mapping shell script

sudo chmod 774 ./Rhapsody_analysis/shell_scripts/Rhapsody_mapping.sh

#Mapping

Go on to the ./Rhapsody_analysis directory on the command line, and 
type ./shell_scripts/Rhapsody_mapping.sh with appropriate options.
mapping result is exported on the ./results directory.

Usage: Rhapsody_mapping.sh [OPTIONS...]
 -h |--help              Display help
 --threads VALUE         Number of CPU threads to use. Default is to use all threads.
 --dataID VALUE          Data ID of read1 and read2 fastq.gz files. read1 (dataID_S[XX]_L00X_R1_001.fastq.gz) is cell barcode read,
                         and read2 (dataID_S[XX]_L00X_R2_001.fastq.gz) is cDNA or tag read.
 --species VALUE         Specify sample species name. Acceptable values are mmu or hsa or rat or macaca
 --library_type VALUE    Specify library type. Acceptable values are WTA, targeted, hashtag, sampletag, ADT, or Abseq.
                         set hashtag if you used BioLegend Hashtags, set sampletag if you used BD Sampletags, set ADT if you use 
                         BioLegend Totalseq antibodies or streptavidin, set Abseq if you used BD Abseq antibodies.
 --expect_cells VALUE    Specify expected cell number of the library. Default=20000. Only uses WTA or targeted mode.
 --index VALUE           Specify index file or folder path. Specify STAR index (for WTA data) or bowtie2 index (hashtag, sampletag, ADT, or Abseq data).
 --BAM VALUE             BAM flag (0=do not save BAM file, 1=save BAM file, default=0)
 --in_dir VALUE          PATH to input .fastq.gz data-stored directory
 --out_dir VALUE         PATH to output directory for saving concatenated raw .fastq.gz data.
 --downsample VALUE      Specify downsampling raw read number if you want to perform downsampling analysis. Default=0 
                         (do not perform downsampling analysis)
 --UMI                   Add this option if you want to perform UMI compression by BD Rhapsody 8-base UMI. Default=OFF (no UMI compression)
