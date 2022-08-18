#!/usr/bin/env Rscript
#extract mapping statistics from SAM files

fnames = dir(pattern = "*_cutadapt_summary.log")
num_of_data = length(fnames)
filename = commandArgs(trailingOnly=TRUE)[1] 

res=read.table(fnames[1], header=TRUE, sep="\t", quote="")

for (i in 2:num_of_data){
 tmp = read.table(fnames[i], header=TRUE, sep="\t", quote="")
 res=rbind(res,tmp)
 }

res=res[,2:ncol(res)]

tmp = colSums(res)
res=rbind(res,tmp)
res=res[nrow(res),c(1,6), drop=F]

invisible(lapply(fnames, file.remove))

fnames = gsub("\\..+$", "", fnames) # remove ".log" from file names
fnames = gsub("_cutadapt_summary", "", fnames) 
colnames(res)=c("total_reads", "passed_reads")

file.name = sprintf("./result/%s_cutadapt.log", filename)
write.table(res, file.name, row.names=T, sep="\t", quote=F)
