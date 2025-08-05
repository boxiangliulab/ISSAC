site_MI<-read.table("/home/users/nus/e0950183/scratch/new_simulation_result/test_pheno/sample400_iso2_cell100_tested_6_27.site")
site_IR<-read.table("/home/users/nus/e0950183/scratch/new_simulation_result/test_pheno/sample400_iso2_cell100_tested_IR_6_27.site")
total_site<-rbind(site_MI,site_IR)
intron<-read.table("/home/users/nus/e0950183/scratch/new_simulation_result/test_pheno/sample400_iso2_cell100_tested_6_27.intron")


genotype<-read.table("/data/projects/11003054/e0950183/compare_site_intron_based/effect_size/genotype/simu_geno.txt",sep=" ",header=TRUE)

sample_num<-c(50,100,200,400)

library(data.table)

ground_truth<-read.table("/data/projects/11003054/e0950183/compare_site_intron_based/effect_size/result/ground_truth_iso2.txt",sep="\t",header=FALSE)
ground_truth$chr<-"NA"
ground_truth$start<-"NA"
ground_truth$end<-"NA"

for(i in 1:nrow(ground_truth)){
    ground_truth$chr[i]<-strsplit(ground_truth[i,3],"_")[[1]][1]
    ground_truth$start[i]<-strsplit(ground_truth[i,3],"_")[[1]][2]
    ground_truth$end[i]<-strsplit(ground_truth[i,3],"_")[[1]][3]
}

for(i in 1:4){
    ###ISSAC pheno
    ISSAC_pheno<-fread(paste("/home/users/nus/e0950183/scratch/new_simulation_result/phenotype/sample400_iso2_effect1_5.prop",sep=""))
    ##ISSAC_pheno<-fread(paste("/home/users/nus/e0950183/scratch/new_simulation_result/metacell/pheno/sample400_iso2_effect1_5_meta2_cell20.prop",sep=""))
    ISSAC_pheno<-as.data.frame(ISSAC_pheno)
    rownames(ISSAC_pheno)<-ISSAC_pheno$V1
    ISSAC_pheno<-ISSAC_pheno[site_MI$V1,]
    ###ISSAC IR pheno
    ISSAC_IR_pheno<-fread(paste("/home/users/nus/e0950183/scratch/new_simulation_result/IR/sample400_iso2_cell100_effect_1_5/sample400_meta2_IR.prop",sep=""))
    ##ISSAC_IR_pheno<-fread(paste("/home/users/nus/e0950183/scratch/new_simulation_result/IR/sample400_iso2_cell20_effect_1_5/sample400_meta2_IR.prop",sep=""))
    ISSAC_IR_pheno<-as.data.frame(ISSAC_IR_pheno)
    rownames(ISSAC_IR_pheno)<-ISSAC_IR_pheno$V1
    ISSAC_IR_pheno<-ISSAC_IR_pheno[site_IR$V1,]
    ISSAC_total<-rbind(ISSAC_pheno,ISSAC_IR_pheno)
    ISSAC_total<-ISSAC_total[,-1]
    sample_ISSAC1<-paste("S",c((1+100):(sample_num[i]+100)),":1",sep="")
    sample_ISSAC2<-paste("S",c((1+100):(sample_num[i]+100)),":2",sep="")
    ISSAC_total<-ISSAC_total[,c(sample_ISSAC1,sample_ISSAC2)]
    ISSAC_result<-as.data.frame(matrix(NA,0,3))
    for(j in 1:nrow(ISSAC_total)){
         chr<-strsplit(rownames(ISSAC_total)[j],":")[[1]][1]
         site<-strsplit(rownames(ISSAC_total)[j],":")[[1]][3]
         num_1<-intersect(which(ground_truth$chr==chr),which(as.numeric(ground_truth$start)<as.numeric(site)))
         num_2<-intersect(which(as.numeric(ground_truth$end)>as.numeric(site)),num_1)
         if((length(num_2)>0)&&(chr!="chrX")&&(chr!="chrY")){
            ground_truth_snp<-unique(ground_truth[num_2,2])
            tmp_pheno<-ISSAC_total[j,]
            tmp_geno<-ISSAC_total[j,]
            for(k in 1:ncol(tmp_geno)){
                sample<-strsplit(colnames(tmp_geno)[k],":")[[1]][1]
                tmp_geno[1,k]<-genotype[ground_truth_snp,sample]
            }
            assoc<-cor.test(as.numeric(tmp_geno[1,]),as.numeric(tmp_pheno[1,]))$p.value
            tmp<-as.data.frame(matrix(NA,1,2))
            tmp[1,1]<-rownames(ISSAC_total)[j]
            tmp[1,2]<-ground_truth_snp
            tmp[1,3]<-assoc
            ISSAC_result<-rbind(ISSAC_result,tmp)
            }}
            
    ###leafcutter pheno
    leafcutter_pheno<-fread(paste("/home/users/nus/e0950183/scratch/new_simulation_result/phenotype/sample400_iso2/sample400_iso2_effect1_5.bed.gz",sep=""),sep="\t")
    ##leafcutter_pheno<-fread(paste("/home/users/nus/e0950183/scratch/new_simulation_result/phenotype/sample400_iso2_effect1_5_cell20.bed.gz",sep=""),sep="\t")
    leafcutter_pheno<-as.data.frame(leafcutter_pheno)
    rownames(leafcutter_pheno)<-leafcutter_pheno$ID
    leafcutter_pheno<-leafcutter_pheno[intron$V1,]
    subsample<-c(paste("S",c((1+100):(sample_num[i]+100)),sep=""))
    ### measure association
    leafcutter_result<-as.data.frame(matrix(NA,0,3))
    for(j in 1:nrow(leafcutter_pheno)){
        chr<-strsplit(rownames(leafcutter_pheno)[j],":")[[1]][1]
        start_intron<-strsplit(rownames(leafcutter_pheno)[j],":")[[1]][3]
        end_intron<-strsplit(rownames(leafcutter_pheno)[j],":")[[1]][3]
        num_1<-intersect(which(ground_truth$chr==chr),which(as.numeric(ground_truth$start)<as.numeric(start_intron)))
        num_2<-intersect(which(as.numeric(ground_truth$end)>as.numeric(end_intron)),num_1)
         if((length(num_2)>0)&&(chr!="chrX")&&(chr!="chrY")){
            ground_truth_snp<-unique(ground_truth[num_2,2])
            tmp_pheno<-leafcutter_pheno[j,subsample]
            tmp_geno<-leafcutter_pheno[j,subsample]
            for(k in 1:ncol(tmp_geno)){
                sample<-colnames(tmp_pheno)[k]
                tmp_geno[1,k]<-genotype[ground_truth_snp,sample]
            }
            assoc<-cor.test(as.numeric(tmp_geno[1,]),as.numeric(tmp_pheno[1,]))$p.value
            tmp<-as.data.frame(matrix(NA,1,3))
            tmp[1,1]<-rownames(leafcutter_pheno)[j]
            tmp[1,2]<-ground_truth_snp
            tmp[1,3]<-assoc
            leafcutter_result<-rbind(leafcutter_result,tmp)
            }
    }
    
    ###Check number of sGenes
    #ISSAC_result$FDR<-p.adjust(ISSAC_result$V3,method="fdr",n=nrow(ISSAC_result))
    #num<-which(ISSAC_result$FDR<0.05)
    #print(length(unique(ISSAC_result$V2[num])))
    #leafcutter_result$FDR<-p.adjust(leafcutter_result$V3,method="fdr",n=nrow(leafcutter_result))
    #num<-which(leafcutter_result$FDR<0.05)
    #print(length(unique(leafcutter_result$V2[num])))


    for(k in 1:nrow(ISSAC_result)){
        ISSAC_result$FDR[k]<-p.adjust(ISSAC_result$V3[k],method="fdr",n=1000)
    }
    for(k in 1:nrow(leafcutter_result)){
        leafcutter_result$FDR[k]<-p.adjust(leafcutter_result$V3[k],method="fdr",n=1000)
    }
    num<-which(ISSAC_result$FDR<0.05)
    print(length(unique(ISSAC_result$V2[num])))
    num<-which(leafcutter_result$FDR<0.05)
    print(length(unique(leafcutter_result$V2[num])))
}


num<-which(ISSAC_result$V3<1e-2)
print(length(unique(ISSAC_result$V2[num])))

num<-which(leafcutter_result$V3<1e-2)
