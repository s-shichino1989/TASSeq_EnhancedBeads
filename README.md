# TAS-Seq mapping pipeline
Mapping pipeline for TAS-Seq / BD Rhapsody WTA data (only for BD Rhapsody enhanced beads)  

To use this pipline, please clone repository, extract, and change folder name as "Rhapsody_analysis".  

Please first see setup.sh to setup analytical environment.  
This package is tested for Ubuntu 20.04 LTS, python3.8, pip, R-4.2.1.  

You could setup analytical environment by moving on to the Rhapsody_analysis directory and type as below.  

```bash
sudo sh setup.sh
```

In the setup.sh script, ask whether re-install R 4.2.1, install required python modules/R packages automatically.  
Please select and type y/Y(Yes) or n/N(No).  

If you do not perform automatic install of modules/packages, you could check required modules/packages  
manually by seeing ./source_file/requirements.txt (for python module installation via pip)  
and ./Rscripts/setup.R (R package installation).  
  
Once setup.sh is finished, add permission to main mapping shell script.  

```bash
sudo chmod 774 ./Rhapsody_analysis/shell_scripts/Rhapsody_mapping.sh
```

## Brief explanation of the pipeline workflow

* Remove adapters, 5' polyG, and 3' polyA streches by Cutadapt 4.1. Quality filtering also performed by Cutadapt 4.1.
* Base composition analysis after quality filtering/trimming by FastQC and Seqkit.
* Perform mapping and counting by STARsolo 2.7.10a (for WTA and targeted reads) or Bowtie2-2.4.5 (for hashtag/sampletag/Totalseq/Abseq reads)
* small-RNA annotations were removed from reference gtf file to build STAR index files.
* Required index files are automatically built by running setup.sh script.
* Identify valid cell barcodes by using DropletUtils and Dropkick packages and export count matrix data of survived cells (.txt.gz file)
* Export spliced- and un-spliced count data for further RNA velocity analysis.
* Export mapping report html file that contains basic statistics of mapping results with associated figures.
* Each cell barcode are exported as 27nt nucleotide bases without two 4-base spacers.

## Dependencies
STARsolo 2.7.10a (https://github.com/alexdobin/STAR)  
Cutadapt 4.1 (https://cutadapt.readthedocs.io/en/stable/)  
Seqkit 2.2.0 (https://bioinf.shenwei.me/seqkit/)  
FastQC v0.11.9 (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
Bowtie2 v2.4.5 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)  
R 4.2.1 (https://cran.r-project.org/) and associated packages (please see ./Rsciprs/setup.R)  
Python3 (3.8 or higher), pip, and associated packages (please see ./source_file/requirements.txt)  

## Usage

Go on to the ./Rhapsody_analysis directory on the command line, and 
type ./shell_scripts/Rhapsody_mapping.sh with appropriate options.  
Mapping result is exported in the ./results directory.

```
Usage: Rhapsody_mapping.sh [OPTIONS...]
 -h |--help                Display help
 
 --threads VALUE           Number of CPU threads to use. Default is to use all threads.
 
 --dataID [required]       Data ID of read1 and read2 fastq.gz files. read1 (dataID_S[XX]_L00X_R1_001.fastq.gz) is cell barcode read,
                           and read2 (dataID_S[XX]_L00X_R2_001.fastq.gz) is cDNA or tag read.
                         
 --species [required]      Specify sample species name. Acceptable values are mmu or hsa or rat or macaca
 
 --library_type [required] Specify library type. Acceptable values are WTA, targeted, hashtag, sampletag, ADT, or Abseq.
                           set hashtag if you used BioLegend Hashtags, set sampletag if you used BD Sampletags, set ADT if you use 
                           BioLegend Totalseq antibodies or streptavidin, set Abseq if you used BD Abseq antibodies.
                           scTCR/scBCR-seq mode is under development.
                         
 --expect_cells VALUE      Specify expected cell number of the library. Default=20000. Only uses WTA or targeted mode.
 
 --index [required]        Specify index file or folder path. Specify STAR index (for WTA data) 
                           or bowtie2 index (hashtag, sampletag, ADT, or Abseq data).
                         
 --BAM VALUE               BAM flag (0=do not save BAM file, 1=save BAM file, default=0)
 
 --in_dir [required]       PATH to input .fastq.gz data-stored directory
 
 --out_dir [required]      PATH to output directory for saving concatenated raw .fastq.gz data.
 
 --downsample VALUE        Specify downsampling raw read number if you want to perform downsampling analysis. 
                           Default=0 (do not perform downsampling analysis)
                         
 --UMI                     Add this option if you want to perform UMI compression by BD Rhapsody 8-base UMI. 
                           Default=OFF (no UMI compression)
```

## Reference
TAS-Seq is a robust and sensitive amplification method for bead-based scRNA-seq.  
Communications Biology volume 5, 602, 2022(https://www.nature.com/articles/s42003-022-03536-0)

## Author
Shigeyuki Shichino (s_shichino@rs.tus.ac.jp)

## Licence
GPL-3
