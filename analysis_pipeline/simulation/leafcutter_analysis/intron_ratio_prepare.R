effect<-"1_7"
intron<-read.table(paste("/home/users/nus/e0950183/scratch/new_simulation_result/phenotype/sample400_iso2/sample400_iso2_effect",effect,".intron.out",sep=""),sep=" ")
splicesite<-read.table(paste("/home/users/nus/e0950183/scratch/new_simulation_result/phenotype/sample400_iso2/sample400_iso2_effect",effect,".site",sep=""),sep=" ",header=TRUE)

newintron<-as.data.frame(matrix(NA,0,ncol(intron)))

for(i in 1:nrow(intron)){
    clu<-intron[i,1]
    num<-which(intron$V1==clu)
    if(length(num)>=2){
        total=rep(0,ncol(intron)-2)
        for(n in 1:length(num)){
            tmp<-as.numeric(intron[num[n],3:ncol(intron)])
            tmp_total=tmp+total
            total=tmp_total
        }
        ratio=as.numeric(intron[i,3:ncol(intron)])/total
        NA_num<-which(is.na(ratio))
        if(length(NA_num)/length(ratio)<0.6){
        denominater<-as.numeric(intron[i,3:ncol(intron)])[-NA_num]
        tmp_middle=round(median(denominater))/round(median(total[-NA_num]))
        ratio[NA_num]<-tmp_middle
        tmp_intron<-as.data.frame(matrix(NA,1,ncol(intron)))
        tmp_intron[1,1:2]<-intron[i,1:2]
        tmp_intron[1,3:ncol(intron)]<-ratio
        newintron<-rbind(newintron,tmp_intron)}
    }
}

colnames(newintron)[1:2]<-c("clu","intron")
colnames(newintron)[3:ncol(newintron)]<-colnames(splicesite)[1:(ncol(splicesite))]

write.table(newintron,paste("/home/users/nus/e0950183/scratch/new_simulation_result/phenotype/sample400_iso2/sample400_iso2_effect",effect,".intron.pheno",sep=""),sep=" ",row.names=FALSE,col.names=TRUE,quote=FALSE)
