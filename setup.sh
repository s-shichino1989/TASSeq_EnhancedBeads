#!/bin/bash -i
#setup script for analysis environment for BD Rhapsody TAS-Seq data
#required software (for Ubuntu20.04).
#If you use CentOS7 substitute apt to yum and adjust package names)


sudo apt-get update -y
sudo apt-get upgrade -y
sudo apt-get install -y parallel wget curl pigz
sudo apt-get install -y libncurses5-dev
sudo apt-get install -y zlib1g-dev
sudo apt-get install -y libbz2-dev
sudo apt-get install -y liblzma-dev
sudo apt-get install -y gcc
sudo apt-get install -y libtool texinfo dpkg-dev pkg-config build-essential
sudo apt-get install -y make
sudo apt-get install -y mpich
sudo apt-get install -y python3-pip libssl-dev libpcre3 libpcre3-dev
sudo apt-get install -y nim libopenblas-base
sudo apt-get install -y hdf5-tools hdf5-helpers libhdf5-dev libhdf5-doc libhdf5-serial-dev
sudo apt-get install -y libcurl4-openssl-dev libxml2-dev
sudo apt-get install -y libgmp3-dev
sudo apt-get install -y gdebi
sudo apt-get install -y libgsl-dev default-jre libgeos++-dev libgeos-dev libgeos-doc
sudo apt-get install -y libhts-dev libboost-all-dev libmagick++-dev libboost-dev
sudo apt-get autoremove -y


echo -n "Do you remove current R and re-install clean R 4.2.1? [y/N]: "
read ANS
 
case $ANS in
  [Yy]* )

    echo "Remove current R if already installed and re-install clean R 4.2.1."
    sudo apt-get purge -y r-base* r-cran-* r-recommended
    sudo apt-get autoremove -y
    sudo rm -R /usr/local/lib/R/site-library
    #install R 4.2.1
    sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
    sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
    sudo apt-get update -y
    sudo apt-get upgrade -y
    sudo apt-get autoremove -y
    sudo apt-get install -y r-base r-base-dev r-recommended
    ;;
  * )
    echo "Skipped re-install of clean R 4.2.1."
    ;;
esac

#installed from source and add PATH to local environment

#STAR
wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz
tar -xzf 2.7.10a.tar.gz
cd STAR-2.7.10a
cd STAR/source
make -j 64 STAR
cd ../../../
sudo cp -r ./STAR-2.7.10a/ /usr/local/STAR-2.7.10a/
sudo ln -fs /usr/local/STAR-2.7.10a/STAR/source/STAR /usr/bin
sudo rm -Rf ./STAR-2.7.10a/

#samtools
if type "samtools" > /dev/null 2>&1
then
 echo "samtools is already installed"
else
 wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
 tar -jxvf samtools-1.9.tar.bz2
 cd ./samtools-1.9/
 ./configure --prefix=/usr/local/
 make -j 8
 sudo make install
 cd ..
 rm samtools-1.9.tar.bz2
 rm -Rf ./samtools-1.9/
fi

#bowtie2 (add symbolic link)
if type "bowtie2" > /dev/null 2>&1
then
 echo "bowtie2 is already installed"
else
 wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.5/bowtie2-2.4.5-linux-x86_64.zip
 unzip bowtie2-2.4.5-linux-x86_64.zip
 sudo mv bowtie2-2.4.5-linux-x86_64 /usr/local/bowtie2-2.4.5
 sudo ln -fs /usr/local/bowtie2-2.4.5/bowtie2 /usr/bin
 sudo ln -fs /usr/local/bowtie2-2.4.5/bowtie2-build /usr/bin
 rm bowtie2-2.4.5-linux-x86_64.zip
fi

#FastQC (add symbolic link)
if type "fastqc" > /dev/null 2>&1
then
 echo "fastqc is already installed"
else
 wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
 unzip fastqc_v0.11.9.zip
 sudo chmod 774 ./FastQC/fastqc
 sudo mv FastQC /usr/local/FastQC
 sudo ln -fs /usr/local/FastQC/fastqc /usr/bin
 rm fastqc_v0.11.9.zip
fi

#seqkit (add symbolic link)
if type "seqkit" > /dev/null 2>&1
then
 echo "seqkit is already installed"
else
 wget https://github.com/shenwei356/seqkit/releases/download/v2.2.0/seqkit_linux_amd64.tar.gz
 tar zxvf seqkit_linux_amd64.tar.gz
 sudo mv seqkit /usr/local/seqkit
 sudo ln -fs /usr/local/seqkit /usr/bin
 rm seqkit_linux_amd64.tar.gz
fi

#fftw3
wget http://www.fftw.org/fftw-3.3.9.tar.gz
tar zxvf fftw-3.3.9.tar.gz
cd fftw-3.3.9
./configure --prefix=/usr/local/fftw3 CC=gcc MPICC=mpicc F77=gfortran --enable-mpi --enable-threads --enable-shared --enable-static
make 
sudo make install
cd ..
rm fftw-3.3.9.tar.gz
rm -R fftw-3.3.9
sudo apt-get install -y libfftw3-dev



#install rstudio
if type "rstudio" > /dev/null 2>&1
then
 echo "rstudio is already installed"
else
 wget https://download1.rstudio.org/desktop/bionic/amd64/rstudio-2022.07.1-554-amd64.deb
 sudo gdebi -n rstudio-2022.07.1-554-amd64.deb
 rm rstudio-2022.07.1-554-amd64.deb
fi

#required python3 modules (pip3)
#python3.8 or higher
echo -n "Do you install required python modules via pip? [y/N]: "
read ANS
 
case $ANS in
  [Yy]* )
    python3 -m pip install --user --upgrade pip
    python3 -m pip install -r --user ./source_file/requirements.txt
    python3 -m pip install --user ./Rhapsody_python_new/
    ;;
  * )
    echo "Skipping install of required python modules. Please see ./source_file/requirements.txt to check required python modules manually."
    ;;
esac

# install pandoc
if type "pandoc" > /dev/null 2>&1
then
 echo "pandoc is already installed"
else
 wget https://github.com/jgm/pandoc/releases/download/2.11.4/pandoc-2.11.4-1-amd64.deb
 sudo gdebi -n pandoc-2.11.4-1-amd64.deb
 rm pandoc-2.11.4-1-amd64.deb
fi


#install required R packages
echo -n "Do you install required R packages automatically? [y/N]: "
read ANS
 
case $ANS in
  [Yy]* )
    Rscript ./Rscripts/setup.R
    ;;
  * )
    echo "Skipping install of required R packages. Please see ./Rscripts/setup.R to check required R packages manually."
    ;;
esac

#find out python3 path and modify source R scripts for appropriate PATH
tmp=`which python3`
cat ./Rscripts/library_source_Seurat_template.R | sed -e 's/path_to_python3/${tmp}/g' > ./Rscripts/library_source_Seurat.R

#create bowtie2 index for tag libraries
sh ./shell_scripts/reference_build_cellHashing.sh

#create STAR index
cd ./reference/fasta_reference_human
sh human_STAR_build.sh
cd ../fasta_reference_mouse
sh mouse_STAR_build.sh
cd ../fasta_reference_rat
sh rat_STAR_build.sh
cd ../fasta_reference_macaca
sh macaca_STAR_build.sh
cd ../../

