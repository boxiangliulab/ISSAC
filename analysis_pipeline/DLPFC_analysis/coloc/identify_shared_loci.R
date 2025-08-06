args <- commandArgs(trailingOnly = TRUE)
cell <- as.numeric(args[1])

test <-as.numeric(args[2])

cell_type<-c("astrocytes", "Excitatory", "Inhibitory", "microglia", "oligodendrocyte", "OPC", "vascular")

summary<-read.table(paste("/home/e0950183/project/brain_sQTL/coloc/sig_splicesite/analysis/",cell_type[cell],".add_gene_sum",sep=""),sep="\t")
num<-which(summary$V7=="sig_splice")
newsummary<-summary[num,]
summary<-newsummary

gwas<-c("selected_AD_Acta_Neuropathologica_Communications_hg38_1e_7.txt.gz",
"selected_ALS_NatComm_hg38_1e_7.txt.gz",        
"selected_ALS_Neuron_hg38_1e_7.txt.gz",              
"selected_LBD_NatGenet_hg38_1e_7.txt.gz",          
"selected_PD_Movement_Disorders_hg38_1e_7.txt.gz",
"selected_AD_NatGenet_2019_hg38_1e_7.txt.gz",                         
"selected_ALS_NatGenet_2016_hg38_1e_7.txt.gz",  
"selected_ANSD_NatGenet_hg38_1e_7.txt.gz",           
"selected_Neuroticism_NatGenet_hg38_1e_7.txt.gz",  
"selected_SCZ_Nature_hg38_1e_7.txt.gz",
"selected_AD_NatGenet_2021_hg38_1e_7.txt.gz",                         
"selected_ALS_NatGenet_2021_hg38_1e_7.txt.gz",  
"selected_Brain_Disease_NatGenet_hg38_1e_7.txt.gz",  
"selected_PD_Lancet_Neurology_hg38_1e_7.txt.gz")

library(data.table)
for(i in test:test){
    coloc<-fread(paste("/home/e1101919/project/brain_sQTL/1e-7_GWAS_coloc/",gwas[i],sep=""),sep="\t",header=TRUE)
    coloc<-as.data.frame(coloc)
    split_list <- lapply(coloc$variant_id, function(x) strsplit(x, ":")[[1]])
    df <- data.frame(do.call(cbind, lapply(split_list, `length<-`, max(sapply(split_list, length)))))
    coloc$chr<-as.character(df[1,])
    coloc$pos<-as.character(df[2,])
    summary$shared<-0
    for(j in 1:nrow(summary)){
        site_chr<-strsplit(summary[j,1],":")[[1]][1]
        site_pos<-as.numeric(strsplit(summary[j,1],":")[[1]][3])
        num_chr<-which(coloc$chr==site_chr)
        num_pos<-intersect(which(as.numeric(coloc$pos)<site_pos+500000),which(as.numeric(coloc$pos)>site_pos-500000))
        if((length(num_chr)>0)&&(length(num_pos)>0)){
            num<-intersect(num_chr,num_pos)
            if(length(num)>0)summary$shared[j]<-1
        }
    }
    num<-which(summary$shared==1)
    if(length(num)>0){
        newsum<-summary[num,]
        write.table(newsum,paste("/home/e0950183/project/brain_sQTL/coloc/shared_loci/",cell_type[cell],gwas[i],".txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
    }
}
