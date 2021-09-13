
#title: "just give the gene name"
#output: html_document

  

setwd("/Users/hoangduy/Documents/Data_Analyse_Ting/PAP_Paper_Lena_TingsProject/")
library(pheatmap)
library(readxl)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(BiocManager)
library(matrixStats)
library(stringr)

#Import Dataset
TPM_value_raw <- read_excel("Lena_2020_microglia_bulkRNA_TPM_value.xlsx")

#remove outliers
TPM_value_raw$MGcKO_Ctrl_18_TPM <- NULL
TPM_value_raw$MGcKO_Ctrl_20_TPM <- NULL


#choose the genes
Genes<-c("Scya1")
as.character(Genes)

Genes_data <-data.frame()
plot_list = list()
for (k in Genes){
#take the data from the database
Genes_data<- subset(TPM_value_raw, TPM_value_raw$Gene_symbol == k)
    if (k%in%TPM_value_raw$Gene_symbol == TRUE) {
Genes_data<-as.data.frame((Genes_data))
Genes_data<- na.omit(Genes_data)
rownames(Genes_data)<-Genes_data$Gene_symbol
Genes_data <- Genes_data[-c(1:4)]
Genes_data <- data.matrix(Genes_data)
#calculate the avg and sd from TPM values
MGcKO         <- mean(Genes_data[k,c(1:5)])
MGcKO_sd      <- sd((Genes_data[k,c(1:5)]))
MGcKO_Ctrl    <- mean(Genes_data[k,c(6:8)])
MGcKO_Ctrl_sd <- sd(Genes_data[k,c(6:8)])
MGcKO_AD      <- mean(Genes_data[k,c(9:13)])
MGcKO_AD_sd   <- sd(Genes_data[k,c(9:13)])
MGcKO_AD_Ctrl <- mean(Genes_data[k,c(14:18)])
MGcKO_AD_Ctrl_sd <- sd(Genes_data[k,c(14:18)])
#put the new values in data.frame
Genes_data <- cbind(Genes_data,MGcKO,MGcKO_Ctrl,MGcKO_AD,MGcKO_AD_Ctrl)
Genes_data <- as.data.frame(Genes_data)
Genes_data <- Genes_data[-c(1:18)]
#turn rows-columms around
Genes_data <- t(Genes_data)
Genes_data <- as.data.frame(Genes_data)
#create a columm named Groups
Genes_data$Groups <- rownames(Genes_data)
Groups <- c(as.character(rownames(Genes_data)))
#create a df with sd
SD <- data.frame()
SD <- data.frame(MGcKO_sd,MGcKO_Ctrl_sd,MGcKO_AD_sd,MGcKO_AD_Ctrl_sd)
SD <- t(SD)
SD <- cbind(SD,Groups)

#merge all together
df_merge <- merge(Genes_data,SD,by="Groups")
names(df_merge)[2] <- paste("GoI")
names(df_merge)[3] <- paste("GoI_sd")
df_merge$GoI <- as.numeric(df_merge$GoI)
df_merge$GoI_sd <- as.numeric(df_merge$GoI_sd)
df_merge[,3] <-  as.numeric(df_merge[,3])

#plot timeeeeeee
Cohorts <- factor(df_merge$Groups, level = c('MGcKO', 'MGcKO_Ctrl', 'MGcKO_AD','MGcKO_AD_Ctrl'))
Title_name <- as.character(k)
options(repr.plot.width = 5, repr.plot.height = 5) 
p = ggplot(data = df_merge, aes(x = Cohorts, y = GoI)) +
  geom_errorbar(aes(
    ymin  = GoI - GoI_sd,
    ymax  = GoI + GoI_sd,
    width = 0.15)) +
  geom_line(group=1) + geom_point() +
  ylab("TPM value")+
  labs(title=Title_name) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 8))
plot_list[[k]] = p
print(p)
    } else {
warning(k)
    }
}
print(p)

for (k in Genes) {
  file_name = paste("Cluster_LenaGenes_Check", k, ".tiff", sep="")
  tiff(file_name)
  print(plot_list[[k]])
  dev.off()
}









