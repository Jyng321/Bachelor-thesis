#Hoang Duy Nguyen
#Bachelor thesis
#2021Jul16
#DESeq2 analysis for Lena microglia bulk RNA-seq data
#all samples
#outlier MGcKO_Ctrl 18 and 20
#analysis under r_backup_env

library(dplyr)
library(DESeq2)
library(ggplot2)

#read in data
count.list<-list.files(path = "./featureCounts on collection 55: Counts",
                       pattern = "tabular")
meta.data<-read.csv("./bulk_metadata.csv")
annot<-read.csv("./mouse_gene_ID_annotation.csv",
                stringsAsFactors = FALSE)
raw<-read.csv("./microglia_bulkRNA_allsamples_raw_counts.csv")

count.list
meta.data
head(annot)
head(raw)
colnames(raw)
rownames(raw)<-raw$X
rownames(meta.data)<-meta.data$sample

#####################
sampleTable <- data.frame(sample = meta.data$sample[-c(8,10)],
                         condition = meta.data$genotype[-c(8,10)])
rownames(sampleTable)<-sampleTable$sample
sampleTable

cts<-raw[,rownames(sampleTable)]
head(cts)

#create combinations
(genotype<-as.character(unique(sampleTable$condition)))

combinations<-combn(genotype,2, simplify = F)
combinations[[2]]<-c("MGcKO_AD","MGcKO")
combinations[[3]]<-c("MGcKO_ADCtrl","MGcKO")
combinations[[4]]<-c("MGcKO_AD","MGcKO_Ctrl")
combinations[[5]]<-c("MGcKO_ADCtrl","MGcKO_Ctrl")
class(combinations[[1]])
combinations[[1]][1]

#extraction test
head(cts[which(sampleTable$condition %in% combinations[[1]])])
sampleTable[which(sampleTable$condition %in% combinations[[1]]),]

#first create empty list to prepare sotre the result
stat<-vector(mode = "list", length = 6)
dds_list<-vector(mode = "list", length = 6)

for(i in 1:6){
    #generate DESeq2 object
    dds<-DESeqDataSetFromMatrix(countData = cts[which(sampleTable$condition %in% combinations[[i]])],
                               colData = sampleTable[which(sampleTable$condition %in% combinations[[i]]),],
                               design = ~condition)
    dds
    
    #decide compare directions, always the second one is the ref
    dds$condition<-relevel(dds$condition, ref = combinations[[i]][2])
    
    #run and save DESeq2 analysis
    dds_result<-DESeq(dds)
    dds_list[[i]]<-dds_result
    saveRDS(dds_result,
            file = paste0("./microglia_bulk_DESeq2_",
                          combinations[[i]][1],"_vs_",combinations[[i]][2],
                          ".rds"))
    
    #extract statistic result
    res<-results(dds_result)
    stat[[i]]<-res
}

stat
dds_list

head(annot)
dim(annot)

result<-vector(mode = "list", length = 6)

#    res<-results(dds_list[[1]])
#    res<-res[order(res$padj, decreasing = F),]
#res<-as.data.frame(res)
##    res
#
#rld<-as.data.frame(counts(dds_list[[1]],normalized=T))
#    for (i in 1:ncol(rld)){
#    colnames(rld)[[i]]<-paste0(colnames(rld)[[i]],".","DESeq2normalised")
#}
##rld
#
#    rld$gene_name<-rownames(rld)
#    
#    res$gene_name<-rownames(res) 
#
#    comparison<-paste0(combinations[[1]][1],"_vs_",combinations[[1]][2],"_")
#    colnames(res)[1:6]<-paste0(comparison,
#                              colnames(res)[1:6])
#
##res
##rld
#
#all<-merge(x = res, y = rld, by = "gene_name")
#
#    adjp_position<-grep("_padj",colnames(all))
#    all<-all[order(all[,adjp_position], decreasing = F),]
#
#
#    all$ensembl_ID<-all$gene_name
#    
#    for(m in 1:nrow(all)){
#    position<-match(all$ensembl_ID[[m]], annot$ensembl_gene_id)
#    all$gene_name[[m]]<-annot$mgi_symbol[position]
#}
#    
#    all<-all[,c(ncol(all),c(1:(ncol(all)-1)))]
#    
#    all

#loop script for organizaing results
for (k in 1:6){
    res<-results(dds_list[[k]])
    res<-res[order(res$padj, decreasing = F),]
    res<-as.data.frame(res)

    rld<-as.data.frame(counts(dds_list[[k]],normalized=T))
    
    for (i in 1:ncol(rld)){
    colnames(rld)[[i]]<-paste0(colnames(rld)[[i]],".","DESeq2normalised")
}
    
    rld$gene_name<-rownames(rld)
    res$gene_name<-rownames(res) 
    
    colnames(res)[1:6]<-paste0(combinations[[k]][1],"_vs_",combinations[[k]][2],"_",
                              colnames(res)[1:6])
    
    all<-merge(x = res, y = rld, by = "gene_name")
    
    adjp_position<-grep("_padj",colnames(all))
    all<-all[order(all[,adjp_position], decreasing = F),]
    
    all$ensembl_ID<-all$gene_name
    
    for(m in 1:nrow(all)){
    position<-match(all$ensembl_ID[[m]], annot$ensembl_gene_id)
    all$gene_name[[m]]<-annot$mgi_symbol[position]
}
    all<-all[,c(ncol(all),c(1:(ncol(all)-1)))]
    
    result[[k]]<-all
}

str(result)
#save
for(n in 1:6){
    write.csv(result[[n]],
              file = paste0("microglia_bulk_DESeq2_rankedresult_",
                            combinations[[n]][1],"_vs_",combinations[[n]][2],
                            ".csv"))
}



