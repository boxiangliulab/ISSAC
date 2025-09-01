### First step: select GWAS loci whose p-value < 1e-7 and exclude other GWAS locis whose distance are within 500kb close to the lead GWAS loci.

### Fourth step: perform coloc analysis to the dataset we prepared and obtain the information of PP0-PP4
args <- commandArgs(trailingOnly = TRUE)
trait_seq <- as.numeric(args[1])

cc_trait=c("selected_AD_Acta_Neuropathologica_Communications_hg38.txt",
"selected_AD_NatGenet_2019_hg38.txt",
"selected_AD_NatGenet_2021_hg38.txt", ##
"selected_ALS_NatComm_hg38.txt",##
"selected_ALS_NatGenet_2016_hg38.txt",##
"selected_ALS_NatGenet_2021_hg38.txt",##
"selected_ALS_Neuron_hg38.txt",##
"selected_ANSD_NatGenet_hg38.txt",
"selected_Brain_Disease_NatGenet_hg38.txt",
"selected_LBD_NatGenet_hg38.txt",
"selected_PD_Lancet_Neurology_hg38.txt",
"selected_SCZ_Nature_hg38.txt"
)
gwas<-c("selected_AD_Acta_Neuropathologica_Communications_hg38_1e_7.txt.gz",
"selected_AD_NatGenet_2019_hg38_1e_7.txt.gz",
"selected_AD_NatGenet_2021_hg38_1e_7.txt.gz", 
"selected_ALS_NatComm_hg38_1e_7.txt.gz",      
"selected_ALS_NatGenet_2016_hg38_1e_7.txt.gz",
"selected_ALS_NatGenet_2021_hg38_1e_7.txt.gz",
"selected_ALS_Neuron_hg38_1e_7.txt.gz", 
"selected_ANSD_NatGenet_hg38_1e_7.txt.gz", 
"selected_Brain_Disease_NatGenet_hg38_1e_7.txt.gz",
"selected_LBD_NatGenet_hg38_1e_7.txt.gz",          
"selected_PD_Lancet_Neurology_hg38_1e_7.txt.gz",
"selected_SCZ_Nature_hg38_1e_7.txt.gz")

gwas_quant<-c("selected_Neuroticism_NatGenet_hg38_1e_7.txt.gz",
        "selected_PD_Movement_Disorders_hg38_1e_7.txt.gz"
)
cc_trait_case<-c(3722,71880,90338,1234,12577,29612,20806,31355,31355,2591,18600,74776)
cc_trait_control<-c(1263,383378,1036225,2850,23475,122656,59804,377103,377103,4027,1400000,101023)

quant_trait=c(
        "selected_Neuroticism_NatGenet_hg38.txt",
        "selected_PD_Movement_Disorders_hg38.txt"
)

quant_num=c(449484,17415)

cell_type<-c("vascular","microglia", "OPC", "oligodendrocyte", "astrocytes", "Inhibitory", "Excitatory")
celltype_sample<-c(1296,1917,2287,3649,2615,4833,6546)

library("coloc")
library(dplyr)
library(data.table)

MAF<-fread("/home/e0950183/project/brain_sQTL/coloc/MAF_total.txt.gz",sep="\t")
MAF<-as.data.frame(MAF)
colnames(MAF)<-c("variant_id","maf")

MAF$variant_id<-toupper(MAF$variant_id)
split_list <- lapply(MAF$variant_id, function(x) strsplit(x, ":")[[1]])
df <- data.frame(do.call(cbind, lapply(split_list, `length<-`, max(sapply(split_list, length)))))
chr<-as.character(df[1,])
pos<-as.character(df[2,])
MAF$variant_id<-paste(chr,pos,sep=":")

for(num_1 in trait_seq:trait_seq){
###read in GWAS dataset
GWAS_data<-paste("/home/e1101919/project/brain_sQTL/processed_GWAS/",cc_trait[num_1],sep="")
print("start read in")
GWAS_pre<-fread(GWAS_data,sep="\t",header=TRUE)
print(paste("read in success",cc_trait[num_1]))
##transform upper character
GWAS_pre$variant_id<-toupper(GWAS_pre$variant_id)
##delete all the ref and alt
split_list <- lapply(GWAS_pre$variant_id, function(x) strsplit(x, ":")[[1]])
df <- data.frame(do.call(cbind, lapply(split_list, `length<-`, max(sapply(split_list, length)))))
chr<-as.character(df[1,])
pos<-as.character(df[2,])
GWAS_pre$variant_id<-paste(chr,pos,sep=":")
###read in intron need to be dealed
trait_type="cc"
for(cell in 1:7){
###coloc
file_name<-paste("/home/e0950183/project/brain_sQTL/coloc/shared_loci/",cell_type[cell],gwas[num_1],".txt",sep="")

if(file.exists(file_name)){
sig<-read.table(file_name,sep="\t")
if(nrow(sig)>0){
        test_rs<-as.data.frame(matrix(NA,nrow(sig),14))
        colnames(test_rs)<-c("variant_id","site_id","gene","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","sQTL_pval","gwas_id","gwas_pval","gwas_name","cell_name")
        test_rs$gwas_name<-gwas[num_1]
        test_rs$cell_name<-cell_type[cell]
        test_rs$gene<-sig$V9
        for(j in 1:nrow(sig)){  ###nrow(sig)
                print(j)
                snp_1mb<-fread(paste("/home/e1101919/project/brain_sQTL/coloc_site_sum/",cell_type[cell],"/",sig[j,1],".result",sep=""),sep="\t",header=FALSE)
                colnames(snp_1mb)<-c("site","variant_id","pval_nominal","beta","se")
                snp_1mb<-as.data.frame(snp_1mb)
                snp_1mb$variant_id<-toupper(snp_1mb$variant_id)
                split_list <- lapply(snp_1mb$variant_id, function(x) strsplit(x, ":")[[1]])
                df <- data.frame(do.call(cbind, lapply(split_list, `length<-`, max(sapply(split_list, length)))))
                chr<-as.character(df[1,])
                pos<-as.numeric(df[2,])
                ###restricted to +-500kb
                site_pos<-as.numeric(strsplit(sig[j,2],":")[[1]][2])
                new_snp_num<-intersect(which(pos>=site_pos-500000),which(pos<=site_pos+500000))
                snp_1mb<-snp_1mb[new_snp_num,]
                pos<-as.character(df[2,])
                snp_1mb$variant_id<-paste(chr[new_snp_num],pos[new_snp_num],sep=":")
                snp_1mb<-snp_1mb[order(snp_1mb$pval_nominal),]
                test_rs[j,1]<-snp_1mb[1,2]
                test_rs[j,2]<-snp_1mb[1,1]
                test_rs[j,3]<-sig[j,8]
                test_rs[j,10]<-snp_1mb[1,3]
                input<-merge(snp_1mb,MAF,by="variant_id")
                input<-merge(input,GWAS_pre,by="variant_id",suffixes=c("_sqtl","_gwas"))
                input<-input[order(input$pvalue),]
                se<-input$se_gwas[1]
                test_rs[j,11]<-input$variant_id[1]
                test_rs[j,12]<-input$pvalue[1]
                input$varbeta<-(input$se_gwas)^2
                input$variant_id<-c(1:nrow(input))
                if(nrow(input)>10){
                        if(trait_type=="cc"){
                                num_zero<-which(input$pval_nominal==0)
                                if(length(num_zero)>0)input$pval_nominal[num_zero]<-input$pval_nominal[num_zero]+1e-100
                                result <- coloc.abf(dataset1=list(pvalues=input$pvalue, type="cc", beta=input$beta_gwas,varbeta=input$varbeta,
                                                                  s=cc_trait_case[num_1]/(cc_trait_case[num_1]+cc_trait_control[num_1]),
                                                                  snp = input$variant_id),
                                                    dataset2=list(pvalues=input$pval_nominal, type="quant", N=celltype_sample[cell],snp=input$variant_id), MAF=input$maf)
                                test_rs[j,4:9]<-t(as.data.frame(result$summary))[1,1:6]
                        }
                        if(trait_type=="quant"){
                                result <- coloc.abf(dataset1=list(pvalues=input$pvalue, type="quant", beta=input$beta_gwas,varbeta=input$varbeta, N=quant_num[num_1],snp = input$variant_id),
                          dataset2=list(pvalues=input$pval_nominal, type="quant", N=celltype_sample[cell],snp=input$variant_id), MAF=input$maf)
                                test_rs[j,4:9]<-t(as.data.frame(result$summary))[1,1:6]}}
        }
  if(nrow(test_rs)>0){write.table(test_rs,paste("/home/e0950183/project/brain_sQTL/coloc/result/NG_Cell/",cell_type[cell],cc_trait[num_1],".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)}
}}
}
}
