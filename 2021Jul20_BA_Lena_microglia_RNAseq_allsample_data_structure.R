#Hoang Duy Nguyen
#Bachelor thesis
#2021Jul16
#test DEG run for Lena microglia bulk RNA-seq data
#running environment bulk env

library(DESeq2)
library(dplyr)
library(ggplot2)

#create a list from the raw-data
count.list<-list.files(path = "/Users/hoangduy/Library/Mobile Documents/com~apple~CloudDocs/Uni/Molmed/Ting's Tag/Data_Analyse_Practice_Ting/featureCounts on collection 55: Counts",
                       pattern = "tabular")
count.list
#create a matrix for the table later
meta.data<-matrix(ncol = 3, nrow = 20)

colnames(meta.data)<-c("file_name","seq_number","genotype")

meta.data<-as.data.frame(meta.data)
meta.data

meta.data$file_name<-count.list

strsplit(count.list[1],"_")

(temp<-strsplit(count.list[1],"_")[[1]][1])

gsub("p885s","",temp)
?gsub


for(i in 1:20){
    temp<-strsplit(meta.data$file_name[i],"_")[[1]][1]
    meta.data$seq_number[i]<-gsub("p885s","",temp)
}

meta.data

meta.data$genotype<-c(rep("MGcKO",5),rep("MGcKO_Ctrl",5),
                      rep("MGcKO_AD",5),rep("MGcKO_ADCtrl",5))

meta.data$sample<-NA

meta.data$sample<-paste0(meta.data$genotype,"_",
                         gsub("MG","",meta.data$seq_number))

meta.data

write.csv(meta.data,"/Users/hoangduy/Library/Mobile Documents/com~apple~CloudDocs/Uni/Molmed/Ting's Tag/Data_Analyse_Practice_Ting.csv")

#combine gene raw counts
readtab<-function(x){
    temp<-read.csv(x, sep = "",stringsAsFactors = FALSE)
    rownames(temp)<-temp$Geneid
    temp<-temp[-1]
    return(temp)
}

setwd("/Users/hoangduy/Library/Mobile Documents/com~apple~CloudDocs/Uni/Molmed/Ting's Tag/Data_Analyse_Practice_Ting/featureCounts on collection 55: Counts")
counts.files<-lapply(count.list, readtab)

str(counts.files)

head(counts.files[[1]])

raw<-Reduce(merge, lapply(counts.files, function(x) data.frame(x, rn = row.names(x))))

head(raw)

rownames(raw)<-raw$rn
raw<-raw[,-1]

head(raw)

colnames(raw)

colnames(raw)==meta.data$file_name

for (i in 1:20){
    position<-grep(colnames(raw)[[i]], meta.data$file_name)
    colnames(raw)[[i]]<-meta.data$sample[position]
}

colnames(raw)

head(raw)
colnames(raw)

raw

#data bank
library("org.Mm.eg.db")

essemble<-rownames(raw)
essemble

convert<-mapIds(org.Mm.eg.db, keys = essemble, keytype = "ENSEMBL", column="SYMBOL")
convert
head(convert)
length(convert)

library(biomaRt)
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
annot<-getBM(c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "strand", "start_position", "end_position","gene_biotype"), 
             mart=ensembl)

head(annot)
dim(annot)
write.csv(annot,"/Users/hoangduy/Library/Mobile Documents/com~apple~CloudDocs/Uni/Molmed/Ting's Tag/Data_Analyse_Practice_Ting/mouse_gene_ID_annotation.csv")

raw$gene_symbol<-NA

for(i in 1:nrow(raw)){
    position<-match(as.character(rownames(raw)[[i]]),annot$ensembl_gene_id)
    raw$gene_symbol[[i]]<-annot$mgi_symbol[position]
}

head(raw)

raw<-raw[,c(21,1:20)]
head(raw)
write.csv(raw,"/Users/hoangduy/Library/Mobile Documents/com~apple~CloudDocs/Uni/Molmed/Ting's Tag/Data_Analyse_Practice_Ting/Lena_2020_microglia_bulkRNA_allsamples_raw_counts.csv")

