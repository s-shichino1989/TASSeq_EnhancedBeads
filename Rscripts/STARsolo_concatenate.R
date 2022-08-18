#!/usr/bin/env Rscript
#concatenate STARsolo results

suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(DropletUtils)))
suppressMessages(suppressWarnings(library(zellkonverter)))

hoge = read.table("dir_list.txt", sep="\t", quote="")

filename = commandArgs(trailingOnly=TRUE)[1]
Nlane = commandArgs(trailingOnly=TRUE)[2]
Nlane = as.numeric(Nlane)

#concatenate GeneFull count table
GeneFull = list()

for (i in c(1:Nlane)){
filename1 = as.character(hoge[i,])
GeneFull[[i]] = Read10X(paste("./", filename1,"/", filename1,"/", filename1, "_GeneFull/raw/", sep=""))
}

for (i in c(1:Nlane)){
  if(i > 1){
    features = union(features, rownames(GeneFull[[i]]))
  } else {
    features = rownames(GeneFull[[1]])
  }
}

features = sort(features)

for (i in c(1:Nlane)){
  colnames(GeneFull[[i]])=gsub("_", "", colnames(GeneFull[[i]]))
}

for (i in c(1:Nlane)){
  if(i > 1){
    cells = union(cells, colnames(GeneFull[[i]]))
  } else {
    cells = colnames(GeneFull[[1]])
  }
}

for(i in c(1:Nlane)){
  features_diff = setdiff(features, rownames(GeneFull[[i]]))
  if(length(features_diff)>1){
    hoge = matrix(nrow = length(features_diff), ncol = ncol(GeneFull[[i]]))
    hoge[is.na(hoge)]=0
    hoge = Matrix(hoge, sparse = T)
    rownames(hoge)=features_diff
    colnames(hoge)=colnames(GeneFull[[i]])
    GeneFull[[i]]=rbind(GeneFUll[[i]], hoge)
  }
  cell_diff = setdiff(cells, colnames(GeneFull[[i]]))
  if(length(cell_diff)>1){
    hoge = matrix(nrow = nrow(GeneFull[[i]]), ncol = length(cell_diff))
    hoge[is.na(hoge)]=0
    hoge = Matrix(hoge, sparse = T)
    rownames(hoge)=rownames(GeneFull[[i]])
    colnames(hoge)=cell_diff
    GeneFull[[i]]=cbind(GeneFull[[i]], hoge)
  }
  #reorder matrix
  GeneFull[[i]] = GeneFull[[i]][features, cells]
}

#combine matrix
for(i in c(1:Nlane)){
  if(i > 1){
  GeneFull1 = GeneFull1 + GeneFull[[i]]
  } else {
  GeneFull1 = GeneFull[[1]]
  }
}

#save sparse matrix
save.folder=paste0("./result/", filename, "_results/processed/")
write10xCounts(path=save.folder, x=GeneFull1, barcodes=colnames(GeneFull1), gene.id = rownames(GeneFull1),
               gene.symbol = rownames(GeneFull1), overwrite=TRUE, type="sparse", version="3")

#save h5ad file for dropkick analysis. pickup top 80000 barcode
GeneFull2 = GeneFull1[,order(Matrix::colSums(GeneFull1, na.rm = T), decreasing = T)]
GeneFull2 = GeneFull2[,c(1:80000)]
GeneFull2 = CreateSeuratObject(GeneFull2)
GeneFull2 = as.SingleCellExperiment(GeneFull2)
save.file=paste0("./result/", filename, "_results/processed/", filename, "_counts_top80000BC.h5ad")
writeH5AD(sce=GeneFull2, file=save.file, compression="gzip")

#concatenate Gene count table
Gene = list()

for (i in c(1:Nlane)){
  filename1 = as.character(hoge[i,])
  Gene[[i]] = Read10X(paste("./", filename1,"/", filename1,"/", filename1, "_Gene/raw/", sep=""))
}

for(i in c(1:Nlane)){
  colnames(Gene[[i]])=gsub("_", "", colnames(Gene[[i]]))
}

for(i in c(1:Nlane)){
  features_diff = setdiff(features, rownames(Gene[[i]]))
  if(length(features_diff)>1){
    hoge = matrix(nrow = length(features_diff), ncol = ncol(Gene[[i]]))
    hoge[is.na(hoge)]=0
    hoge = Matrix(hoge, sparse = T)
    rownames(hoge)=features_diff
    colnames(hoge)=colnames(Gene[[i]])
    Gene[[i]]=rbind(Gene[[i]], hoge)
  }
  cell_diff = setdiff(cells, colnames(Gene[[i]]))
  if(length(cell_diff)>1){
    hoge = matrix(nrow = nrow(Gene[[i]]), ncol = length(cell_diff))
    hoge[is.na(hoge)]=0
    hoge = Matrix(hoge, sparse = T)
    rownames(hoge)=rownames(Gene[[i]])
    colnames(hoge)=cell_diff
    Gene[[i]]=cbind(Gene[[i]], hoge)
  }
  #reorder matrix
  Gene[[i]] = Gene[[i]][features, cells]
}

#combine matrix
for(i in c(1:Nlane)){
  if(i > 1){
    Gene1 = Gene1 + Gene[[i]]
  } else {
    Gene1 = Gene[[1]]
  }
}
#save sparse matrix
save.folder=paste0("./result/", filename, "_results/Gene/")
write10xCounts(path=save.folder, x=Gene1, barcodes=colnames(Gene1), gene.id = rownames(Gene1),
               gene.symbol = rownames(Gene1), overwrite=TRUE, type="sparse", version="3")


#concatenate STARsolo summary file
STARsolo_log = list()
for (i in c(1:Nlane)){
  filename1 = as.character(hoge[i,])
  tmp =   read.csv(paste0("./result/", filename, "_results/stats/", filename1, "_Summary.csv"),
                                 header=F, quote="", stringsAsFactors = F)
  if(i>1){
  STARsolo_log_res = cbind(STARsolo_log_res, tmp[,2])
  } else {
  STARsolo_log_res =  tmp[,2,drop=F]
  }
  rownames(STARsolo_log_res) = tmp[,1]
}

STARsolo_log_sum = rowSums(STARsolo_log_res)
STARsolo_log_ave = rowMeans(STARsolo_log_res)

STARsolo_log = cbind(tmp[,1], c(STARsolo_log_sum[1], STARsolo_log_ave[2:length(STARsolo_log_ave)]))
write.csv(STARsolo_log, file = paste0("./result/", filename, "_results/stats/STARsolo_Summary.csv"),
          row.names = F, quote = F)

for (i in c(1:Nlane)){
filename1 = as.character(hoge[i,])
invisible(file.remove(paste0("./result/", filename, "_results/stats/", filename1, "_Summary.csv")))
}
