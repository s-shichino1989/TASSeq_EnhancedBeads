
library(data.table)
library(stringr)

CS1 = read.table("/datadrive/Rhapsody_analysis/reference/BD_CB/BD_CLS1.txt")
CS2 = read.table("/datadrive/Rhapsody_analysis/reference/BD_CB/BD_CLS2.txt")
CS3 = read.table("/datadrive/Rhapsody_analysis/reference/BD_CB/BD_CLS3.txt")

CS1$hash = rownames(CS1)
CS1 = rbind(CS1, c("x", "x"))
CS2$hash = rownames(CS2)
CS2 = rbind(CS2, c("x", "x"))
CS3$hash = rownames(CS3)
CS3 = rbind(CS3, c("x", "x"))

#read annotated R1 csv file
filename = commandArgs(trailingOnly=TRUE)[1] 
hoge = fread(filename, sep=",", header=F)
hoge1 = str_split_fixed(hoge$V1, pattern="-", n=3)
hoge = cbind(hoge1, hoge[,c(2,3)])

fuga = hoge[,1]
colnames(fuga)="hash"
fuga = merge.data.table(fuga, CS1, by="hash", sort = FALSE, all.x = TRUE)

fuga2 = hoge[,2]
colnames(fuga2)="hash"
fuga2 = merge.data.table(fuga2, CS2, by="hash", sort = FALSE, all.x = TRUE)

fuga3 = hoge[,3]
colnames(fuga3)="hash"
fuga3 = merge.data.table(fuga3, CS3, by="hash", sort = FALSE, all.x =TRUE)

hoge = as.data.frame(hoge)

res = data.frame(CLS1=fuga$V1, CLS2=fuga2$V1, CLS3=fuga3$V1, UMI=hoge[,5])
colnames(res)=c("CLS1", "CLS2", "CLS3", "UMI")

#retrive cell barcode and UMI sequence
res1 = paste0(res$CLS1,  res$CLS2, res$CLS3)
res1 = as.data.frame(res1)
res1 = cbind(res1, editops=hoge[,4], UMI=hoge[,5])

#create virtual read quality table
res2 = data.frame(CB_qual=gsub(pattern = ".", replacement = "F", res1[,1]),
                  MR_qual=gsub(pattern = ".", replacement = "F", res1[,3]))
res1 = cbind(res1, res2)

fwrite(res1, file=filename, sep=",", 
       quote=F, row.names=F,col.names=F)

