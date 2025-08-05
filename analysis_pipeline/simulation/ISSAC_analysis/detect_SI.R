iso<-2
sample<-400
effect<-"1_5"
cell<-"100"

intron<-read.table(paste("/home/users/nus/e0950183/scratch/new_simulation_result/IR/sample",sample,"_iso",iso,"_cell",cell,"_effect_",effect,"/single_intron_cluster",sep=""),sep="\t")
total<-read.table(paste("/home/users/nus/e0950183/scratch/new_simulation_result/metacell/pheno/sample",sample,"_iso",iso,"_effect",effect,"_meta2.intron.out",sep=""),sep=" ")
#total<-read.table("/home/users/nus/e0950183/scratch/new_simulation_result/metacell/pheno/sample400_iso2_effect1_5_meta2_cell20.intron.out",sep=" ")
site<-as.data.frame(matrix(NA,0,1))
for(i in 1:nrow(intron)){
    num<-which(total[,1]==intron[i,1])
    chr<-strsplit(total[num,2],":")[[1]][1]
    strand<-strsplit(total[num,2],":")[[1]][2]
    start<-strsplit(total[num,2],":")[[1]][3]
    end<-strsplit(total[num,2],":")[[1]][4]
    site1<-paste(chr,strand,start,sep=":")
    site2<-paste(chr,strand,end,sep=":")
    tmp<-as.data.frame(matrix(NA,2,1))
    tmp[1,1]<-site1
    tmp[2,1]<-site2
    site<-rbind(site,tmp)
}

write.table(site,paste("/home/users/nus/e0950183/scratch/new_simulation_result/IR/sample",sample,"_iso",iso,"_cell",cell,"_effect_",effect,"/single_intron_site",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
