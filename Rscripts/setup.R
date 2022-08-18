#!/usr/bin/env Rscript
#setup R environment
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs = c("data.table", "dplyr", "ggplot2", 
"stringr", "tibble", "tidyr", "mcgv", "voxel", "ggplotify","foreach", "doParallel", "withr",
"Matrix", "mclust", "recommenderlab", "stringdist", "MASS", "future.apply", "future", "doFuture", 
"reticulate", "viridis", "devtools", "factoextra", "flashClust", "cowplot", "XML", "RCurl", "Seurat", 
"SingleCellExperiment","SummarizedExperiment", "GenomicRanges", "org.Mm.eg.db", "org.Hs.eg.db",
 "TCC", "impute", "preprocessCore", "GO.db", "AnnotationDbi", "TCC", "BiocParallel", 
"AnnotationDbi", "biomaRt", "flowDensity", "zellkonverter", "WGCNA", "parallelDist", "ggsci", "filesstrings",
"DropletUtils","flowTrans", "slingshot", "clusterExperiment", "rmdformats", "DT", "ggpubr",
"doRNG", "qs", "dbplyr", "magick", "circlize", "ggnewscale", "msigdbr", "tidyverse", 
"minidown", "extrafont",
"monocle", "multtest", "clusterProfiler", "pcaMethods", "BiocNeighbors", "SingleR", "ComplexHeatmap")

BiocManager::install(pkgs)

#if you use CentOS7 you must install hdf5 manually
#wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.20/src/hdf5-1.8.20.tar.gz
#unpigz hdf5-1.8.20.tar.gz
#tar xf hdf5-1.8.20.tar.gz
#cd hdf5-1.8.20
#./configure --enable-fortran --prefix=/usr/local/hdf5 --enable-cxx --with-pthread=/usr/include/,/usr/lib/x86_64-linux-gnu
#make -j 24
#make check
#sudo make install
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/hdf5/lib
#install.packages("hdf5r", configure.args="--with-hdf5=/usr/local/hdf5/bin/h5cc")
#install.packages("./source_file/SDMTools_1.1-221.2.tar.gz", type="source", repos=NULL, lib="/usr/local/lib/R/site-library")


library(devtools)
install_github('barkasn/fastSave')
install.packages("./source_file/tradeSeq-master.tar.gz", type="source", repos=NULL)
install.packages("./source_file/rDBEC_0.3.2.tar.gz", type="source", repos=NULL)
install_github("velocyto-team/velocyto.R")

