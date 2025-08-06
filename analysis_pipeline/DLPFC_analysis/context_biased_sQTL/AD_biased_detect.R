### read in phenotype
library(lme4)
library(QuantPsyc)
cell<-"vascular"

dat<-read.table(paste("/home/users/nus/e0950183/scratch/brain_major_sub/result/analysis/",cell,".add_gene_sum",sep=""),sep="\t")

num<-which(dat[,7]=="sig_splice")

newdat<-dat[num,]

pheno_name<-read.table(paste("/data/projects/11003054/e0950183/brain_sQTL_AD_sex_cell_biased/AD_biased/",cell,"_pheno_name.txt",sep=""))
geno_name<-read.table("/data/projects/11003054/e0950183/brain_sQTL_AD_sex_cell_biased/AD_biased/geno_name.txt",sep="\t")

PC<-read.table(paste("/home/users/nus/e0950183/scratch/brain_major_sub/PC/",cell,".PC",sep=""),sep="\t",header=TRUE,check.names=FALSE)
final_result<-as.data.frame(matrix(NA,nrow(newdat),11))
final_result$V1<-newdat$V1
final_result$V2<-newdat$V2
final_result$V3<-newdat$V9
final_result$V4<-newdat$V3

for(i in 1:nrow(newdat)){
    print(i)
    pheno<-read.table(paste("/data/projects/11003054/e0950183/brain_sQTL_AD_sex_cell_biased/AD_biased/phenotype/",cell,"/",newdat[i,1],".txt",sep=""),sep=" ",row.names=1)
    colnames(pheno)<-pheno_name$V1
    geno<-read.table(paste("/data/projects/11003054/e0950183/brain_sQTL_AD_sex_cell_biased/AD_biased/genotype/",cell,"/",newdat[i,2],".txt",sep=""),sep="\t")
    pos<-strsplit(newdat[i,2],":")[[1]][2]
    num_pos<-which(geno[,2]==pos)
    geno<-geno[num_pos,]
    colnames(geno)<-geno_name[1,]
    ##split pheno & geno
    total<-as.data.frame(matrix(NA,4,ncol(pheno)))
    for(j in 1:ncol(pheno)){
        total[1,j]<-strsplit(pheno[1,j],":")[[1]][1]
        total[2,j]<-strsplit(pheno[1,j],":")[[1]][2]
        total[3,j]<-strsplit(colnames(pheno)[j],":")[[1]][1]
        tmp_geno<-geno[1,total[3,j]]
        if(startsWith(tmp_geno,'0|0')){total[4,j]<-0}
        if(startsWith(tmp_geno,'0|1')){total[4,j]<-1}
        if(startsWith(tmp_geno,'1|0')){total[4,j]<-1}
        if(startsWith(tmp_geno,'1|1')){total[4,j]<-2}
    }
    colnames(total)<-pheno_name$V1
    ###reformat results
    model_data<-data.frame(suc=as.numeric(total[1,]),
    tot=as.numeric(total[2,]),
    sPC1=as.numeric(PC[1,]),
    sPC2=as.numeric(PC[2,]),
    sPC3=as.numeric(PC[3,]),
    sPC4=as.numeric(PC[4,]),
    sPC5=as.numeric(PC[5,]),
    sPC6=as.numeric(PC[6,]),
    sPC7=as.numeric(PC[7,]),
    sPC8=as.numeric(PC[8,]),
    gPC1=as.numeric(PC[9,]),
    gPC2=as.numeric(PC[10,]),
    gPC3=as.numeric(PC[11,]),
    gPC4=as.numeric(PC[12,]),
    gPC5=as.numeric(PC[13,]),
    sex=as.numeric(PC[14,]),
    age=as.numeric(PC[15,]),
    AD=as.numeric(PC[16,]),
    ROS_MAP=as.numeric(PC[17,]),
    PMI=as.numeric(PC[18,]),
    num_cells=as.numeric(PC[19,]),
    educ=as.numeric(PC[20,]),
    data_source=as.numeric(PC[21,]),
    group=factor(as.character(total[3,])),
    geno=as.numeric(total[4,]),
    geno_sex=as.numeric(total[4,])*as.numeric(PC[14,])
    )
    model_data$fal<-model_data$tot-model_data$suc
    model <- glmer(cbind(suc, fal) ~ sPC1 + sPC2 + sPC3 + sPC4 + sPC5 + sPC6 + sPC7 + sPC8 + gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
        sex + age + AD + ROS_MAP + PMI + num_cells + educ + data_source + geno + geno_sex + (1|group), 
               data = model_data, 
               family = binomial(link = "logit"))
    final_result$V5[i]<-summary(model)$coefficients[,'Pr(>|z|)']['geno']
    final_result$V6[i]<-summary(model)$coefficients[,'Pr(>|z|)']['geno_sex']
    num<-which(model_data$sex==1)
    newmodel<-model_data[num,]
    model <- glmer(cbind(suc, fal) ~ sPC1 + sPC2 + sPC3 + sPC4 + sPC5 + sPC6 + sPC7 + sPC8 + gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
        AD + age + ROS_MAP + PMI + num_cells + educ + data_source + geno+ (1|group), 
               data = newmodel, 
               family = binomial(link = "logit"))
    final_result$V8[i]<-summary(model)$coefficients[,'Pr(>|z|)']['geno']
    final_result$V10[i]<-lm.beta(model)['geno']
    num<-which(model_data$sex==2)
    newmodel<-model_data[num,]
    model <- glm(cbind(suc, fal) ~ sPC1 + sPC2 + sPC3 + sPC4 + sPC5 + sPC6 + sPC7 + sPC8 + gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
        AD + age + ROS_MAP + PMI + num_cells + educ + data_source + geno+ (1|group), 
               data = newmodel, 
               family = binomial(link = "logit"))
    final_result$V9[i]<-summary(model)$coefficients[,'Pr(>|z|)']['geno']
    final_result$V11[i]<-lm.beta(model)['geno']
}

final_result$V7<-p.adjust(final_result$V6,method="fdr",n=nrow(final_result))

write.table(final_result,paste("/data/projects/11003054/e0950183/brain_sQTL_AD_sex_cell_biased/sex_biased/",cell,"_result",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
