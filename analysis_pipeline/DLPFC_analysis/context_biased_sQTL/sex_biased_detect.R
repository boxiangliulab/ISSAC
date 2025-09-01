args <- commandArgs(trailingOnly = TRUE)
cell <- as.character(args[1])

test <-as.numeric(args[2])
### read in phenotype
library(lme4)
library(QuantPsyc)
computeRegressionSlope <- function(x, y) {
  # Ensure the vectors are of the same length
  if (length(x) != length(y)) {
    stop("Vectors x and y must have the same size.")
  }
  
  # Calculate the mean of x and y
  x_mean <- mean(x)
  y_mean <- mean(y)
  
  # Calculate the standard deviation of x and y
  x_stddev <- sqrt(mean((x - x_mean)^2))
  y_stddev <- sqrt(mean((y - y_mean)^2))
  
  # Normalize x and y
  x_normalized <- (x - x_mean) / x_stddev
  y_normalized <- (y - y_mean) / y_stddev
  
  # Compute the numerator and denominator for the slope
  numerator <- sum(x_normalized * y_normalized)
  denominator <- sum(x_normalized^2)
  
  # Calculate and return the slope
  return(numerator / denominator)
}

dat<-read.table(paste("/home/users/nus/e0950183/scratch/brain_major_sub/result/analysis/",cell,".add_gene_sum",sep=""),sep="\t")

num<-which(dat[,7]=="sig_splice")

newdat<-dat[num,]

pheno_name<-read.table(paste("/data/projects/11003054/e0950183/brain_sQTL_AD_sex_cell_biased/AD_biased/",cell,"_pheno_name.txt",sep=""))
geno_name<-read.table("/data/projects/11003054/e0950183/brain_sQTL_AD_sex_cell_biased/AD_biased/geno_name.txt",sep="\t")

PC<-read.table(paste("/data/projects/11003054/e0950183/brain_sQTL/PC/PC/",cell,".PC",sep=""),sep="\t",header=TRUE,check.names=FALSE)
final_result<-as.data.frame(matrix(NA,nrow(newdat),11))
final_result$V1<-newdat$V1
final_result$V2<-newdat$V2
final_result$V3<-newdat$V8
final_result$V4<-newdat$V3

start=1000*test-999
end=1000*test
if(end>nrow(newdat)){end=nrow(newdat)}

for(i in start:end){
    print(i)
    pheno<-read.table(paste("/home/users/nus/e0950183/scratch/brain_AD_sex/phenotype/",cell,"/",newdat[i,1],".txt",sep=""),sep=" ",row.names=1)
    colnames(pheno)<-pheno_name$V1
    
    pheno_beta_compute<-read.table(paste("/home/users/nus/e0950183/scratch/brain_major_sub/result/major/",cell,"_model/",newdat[i,1],".middle",sep=""),skip=1,row.names=1)
    pheno_revise<-(pheno_beta_compute[4,]/pheno_beta_compute[3,]) - pheno_beta_compute[2,]

    geno<-read.table(paste("/home/users/nus/e0950183/scratch/brain_AD_sex/genotype/",cell,"/",newdat[i,2],".txt",sep=""),sep="\t")
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
    ###Perform subsample to AD
    AD_male<-intersect(which(model_data$sex==1),which(model_data$AD==1))
    nonAD_male<-intersect(which(model_data$sex==1),which(model_data$AD==0))
    if(length(AD_male)>length(nonAD_male)){
        newAD_male<-AD_male[sample(1:length(AD_male),length(nonAD_male))]
        new_sample_male<-c(newAD_male,nonAD_male)
    }else{
        newnonAD_male<-nonAD_male[sample(1:length(nonAD_male),length(AD_male))]
        new_sample_male<-c(newnonAD_male,AD_male)
    }
    AD_female<-intersect(which(model_data$sex==2),which(model_data$AD==1))
    nonAD_female<-intersect(which(model_data$sex==2),which(model_data$AD==0))
    
    if(length(AD_female)>length(nonAD_female)){
        newAD_female<-AD_female[sample(1:length(AD_female),length(nonAD_female))]
        new_sample_female<-c(newAD_female,nonAD_female)
    }else{
        newnonAD_female<-nonAD_female[sample(1:length(nonAD_female),length(AD_female))]
        new_sample_female<-c(newnonAD_female,AD_female)
    }
    new_sample<-c(new_sample_male,new_sample_female)
    model_data<-model_data[new_sample,]
    pheno_revise<-pheno_revise[new_sample]
    model <- glmer(cbind(suc, fal) ~ sPC1 + sPC2 + sPC3 + sPC4 + sPC5 + sPC6 + sPC7 + sPC8 + gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
        sex + age + AD + ROS_MAP + PMI + num_cells + educ + data_source + geno + geno_sex+(1|group), 
               data = model_data, 
               family = binomial(link = "logit"))
    final_result$V5[i]<-summary(model)$coefficients[,'Pr(>|z|)']['geno']
    final_result$V6[i]<-summary(model)$coefficients[,'Pr(>|z|)']['geno_sex']
    num<-which(model_data$sex==1)
    newmodel<-model_data[num,]
    model <- glmer(cbind(suc, fal) ~ sPC1 + sPC2 + sPC3 + sPC4 + sPC5 + sPC6 + sPC7 + sPC8 + gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
        AD + age + ROS_MAP + PMI + num_cells + educ + data_source + geno+(1|group), 
               data = newmodel, 
               family = binomial(link = "logit"))
    final_result$V8[i]<-summary(model)$coefficients[,'Pr(>|z|)']['geno']
    final_result$V10[i]<-computeRegressionSlope(as.numeric(pheno_revise[num]),as.numeric(newmodel$geno))
    num<-which(model_data$sex==2)
    newmodel<-model_data[num,]
    model <- glmer(cbind(suc, fal) ~ sPC1 + sPC2 + sPC3 + sPC4 + sPC5 + sPC6 + sPC7 + sPC8 + gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + 
        AD + age + ROS_MAP + PMI + num_cells + educ + data_source + geno+(1|group), 
               data = newmodel, 
               family = binomial(link = "logit"))
    final_result$V9[i]<-summary(model)$coefficients[,'Pr(>|z|)']['geno']
    final_result$V11[i]<-computeRegressionSlope(as.numeric(pheno_revise[num]),as.numeric(newmodel$geno))
}

final_result<-final_result[start:end,]
final_result$V7<-p.adjust(final_result$V6,method="fdr",n=nrow(final_result))

write.table(final_result,paste("/home/users/nus/e0950183/scratch/brain_AD_sex/result/",cell,"_result_sex",test,sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
