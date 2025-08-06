args <- commandArgs(trailingOnly = TRUE)
trait_seq <- as.numeric(args[1])

#total<-read.table("/home/e0950183/project/brain_sQTL/coloc/summary_coloc_H4_over_0_75.txt",sep="\t")
total<-read.table("/home/e0950183/project/brain_sQTL/coloc/summary_coloc_2_24.txt",sep="\t")
num<-which(total$V2=="chr17:+:45974384")
total<-total[num,]
gwas_short<-c("selected_AD_Acta_Neuropathologica_Communications_hg38_1e_7.txt.gz",
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
"selected_SCZ_Nature_hg38_1e_7.txt.gz",
"selected_Neuroticism_NatGenet_hg38_1e_7.txt.gz",
"selected_PD_Movement_Disorders_hg38_1e_7.txt.gz")

gwas_list=c("selected_AD_Acta_Neuropathologica_Communications_hg38.txt",
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
"selected_SCZ_Nature_hg38.txt",
"selected_Neuroticism_NatGenet_hg38.txt",
"selected_PD_Movement_Disorders_hg38.txt")

library(data.table)

cc_trait_case<-c(3722,71880,90338,1234,12577,29612,20806,31355,31355,2591,18600,74776)
cc_trait_control<-c(1263,383378,1036225,2850,23475,122656,59804,377103,377103,4027,1400000,101023)
cc_trait<-cc_trait_case+cc_trait_control
quant_num=c(449484,17415)
gwas_num<-c(cc_trait,quant_num)

for(i in trait_seq:trait_seq){
    num<-which(total$V13==gwas_short[i])
    tmp_total<-total[num,]
    if(length(num)>0){
        gwas_file<-fread(paste("/home/e1101919/project/brain_sQTL/processed_GWAS/",gwas_list[i],sep=""),sep="\t")
        gwas_file<-as.data.frame(gwas_file)
        gwas_file$variant_id<-toupper(gwas_file$variant_id)
##delete all the ref and alt
        split_list <- lapply(gwas_file$variant_id, function(x) strsplit(x, ":")[[1]])
        df <- data.frame(do.call(cbind, lapply(split_list, `length<-`, max(sapply(split_list, length)))))
        chr<-as.character(df[1,])
        pos<-as.character(df[2,])
        gwas_file$variant_id<-paste(chr,pos,sep=":")
        for(j in 1:nrow(tmp_total)){
          snp_1mb<-fread(paste("/home/e1101919/project/brain_sQTL/coloc_site_sum/",tmp_total[j,14],"/",tmp_total[j,2],".result",sep=""),sep="\t",header=FALSE)
          colnames(snp_1mb)<-c("site","variant_id","pval_nominal","beta","se")
          snp_1mb<-as.data.frame(snp_1mb)
          snp_1mb$variant_id<-toupper(snp_1mb$variant_id)
          split_list <- lapply(snp_1mb$variant_id, function(x) strsplit(x, ":")[[1]])
          df <- data.frame(do.call(cbind, lapply(split_list, `length<-`, max(sapply(split_list, length)))))
          chr<-as.character(df[1,])
          pos<-as.numeric(df[2,])
          snp_1mb$variant_id<-paste(chr,pos,sep=":")
          input<-merge(snp_1mb,gwas_file,by="variant_id",suffixes=c("_sqtl","_gwas"))
          input$distance<-0
          site_chr<-strsplit(tmp_total[j,2],":")[[1]][1]
          site_pos<-strsplit(tmp_total[j,2],":")[[1]][3]
          vcf<-read.table(paste("/home/e0950183/project/brain_sQTL/coloc/smr/vcf/",site_chr,"_",site_pos,".vcf",sep=""),sep="\t")
          input$SNP<-NA
          input$chr<-NA
          input$BP<-NA
          input$A1<-NA
          input$A2<-NA
          input$Freq<-NA
          for(k in 1:nrow(input)){
            pos<-strsplit(input$variant_id[k],":")[[1]][2]
            num_snp<-which(vcf$V2==pos)
            input$SNP[k]<-vcf$V3[num_snp]
            input$chr[k]<-strsplit(vcf$V1[num_snp],"r")[[1]][2]
            input$BP[k]<-vcf$V2[num_snp]
            input$A1[k]<-vcf$V4[num_snp]
            input$A2[k]<-vcf$V5[num_snp]
            MAF<-strsplit(vcf$V8[num_snp],";AF=")[[1]][2]
            input$Freq[k]<-1-as.numeric(strsplit(MAF,";")[[1]][1])
          }
          input$Probe_bp<-site_pos
          input$Orientation<-strsplit(tmp_total[j,2],":")[[1]][2]
          input$Gene<-input$site
          tmp_gwas<-input[,c(10,13,14,15,6,7,8)]
          tmp_gwas$N<-gwas_num[i]
          tmp_sQTL<-input[,c(10,11,12,13,14,15,2,11,16,18,17,4,5,3)]
          colnames(tmp_gwas)<-c("SNP", "A1","A2","Freq","b","se","p","N")
          colnames(tmp_sQTL)<-c("SNP", "Chr", "BP","A1","A2","Freq","Probe","Probe_Chr","Probe_bp","Gene","Orientation","b","se","p")
          tmp_gwas_dedup<-tmp_gwas[!duplicated(tmp_gwas$SNP),]
          tmp_sQTL_dedup<-tmp_sQTL[!duplicated(tmp_sQTL$SNP),]
          write.table(tmp_gwas_dedup,paste("/home/e0950183/project/brain_sQTL/coloc/smr/smr_GWAS/",tmp_total[j,14],"-",tmp_total[j,2],"-",tmp_total[j,13],sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
          write.table(tmp_sQTL_dedup,paste("/home/e0950183/project/brain_sQTL/coloc/smr/smr_sQTL/",tmp_total[j,14],"-",tmp_total[j,2],"-",tmp_total[j,13],sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
        }
    }
}

###read in GWAS


###read in sQTL
