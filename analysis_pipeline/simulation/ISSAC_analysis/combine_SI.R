iso<-2

sample<-50

cell<-100
effect<-"1_5"
splice_site<-read.table(paste("/home/users/nus/e0950183/scratch/new_simulation_result/IR/sample400_iso2_cell100_effect_",effect,"/single_intron_site",sep=""))
newdat<-as.data.frame(matrix(0,nrow(splice_site),2*sample))
rownames(newdat)<-splice_site$V1

intron<-read.table(paste("/home/users/nus/e0950183/scratch/new_simulation_result/metacell/pheno/sample400_iso2_effect",effect,"_meta2.intron.out",sep=""))
sample_name<-read.table(paste("/home/users/nus/e0950183/scratch/new_simulation_result/sample_list/sample",sample,"_meta.list",sep=""))
colnames(intron)[3:(2+2*sample)]<-sample_name$V1

colnames(newdat)<-sample_name$V1

meta1<-read.table("/home/users/nus/e0950183/scratch/new_simulation_result/metacell/meta2_1.list")
meta2<-read.table("/home/users/nus/e0950183/scratch/new_simulation_result/metacell/meta2_2.list")
for(i in 1:sample){
    print(i)
    nonsplit<-read.table(paste("/home/users/nus/e0950183/scratch/new_simulation_result/IR/sample400_iso2_cell",cell,"_effect_",effect,"/S",i,".nonsplit",sep=""))
    for(j in 1:nrow(nonsplit)){
        site_chr<-strsplit(nonsplit[j,1],":")[[1]][1]
        site_strand<-strsplit(nonsplit[j,1],":")[[1]][3]
        site_pos<-strsplit(nonsplit[j,1],":")[[1]][2]
        site<-paste(site_chr,site_strand,site_pos,sep=":")
        barcode<-strsplit(nonsplit[j,1],":")[[1]][4]
        if(barcode%in%(meta1$V1)){
        num<-which(rownames(newdat)==site)
        name<-paste("S",i,":1",sep="")
        newdat[num,name]<-newdat[num,name]+1}
        if(barcode%in%(meta2$V1)){
        num<-which(rownames(newdat)==site)
        name<-paste("S",i,":2",sep="")
        newdat[num,name]<-newdat[num,name]+1
        }
    }
}

for(i in 1:nrow(newdat)){
    if(sum(as.numeric(newdat[i,]))==0){
        newdat[i,1]<-NA
    }
}
newdat<-na.omit(newdat)

write.table(newdat,paste("/home/users/nus/e0950183/scratch/new_simulation_result/IR/sample400_iso2_cell",cell,"_effect_",effect,"/sample",sample,"_meta2.nonsplit",sep=""),row.names=TRUE,col.names=TRUE,quote=FALSE)

newdat<-read.table(paste("/home/users/nus/e0950183/scratch/new_simulation_result/IR/sample400_iso2_cell",cell,"_effect_",effect,"/sample",sample,"_meta2.nonsplit",sep=""),header=TRUE,check.names=FALSE)
for(i in 1:nrow(intron)){
    chr<-strsplit(intron[i,2],":")[[1]][1]
    strand<-strsplit(intron[i,2],":")[[1]][2]
    start<-strsplit(intron[i,2],":")[[1]][3]
    end<-strsplit(intron[i,2],":")[[1]][4]
    site1<-paste(chr,strand,start,sep=":")
    site2<-paste(chr,strand,end,sep=":")
    intron$start[i]<-site1
    intron$end[i]<-site2
}

###rearrange intron position


newintron<-intron[,colnames(newdat)]
for(i in 1:nrow(newdat)){
    num_1<-which((intron$start==rownames(newdat)[i]))
    num_2<-which((intron$end==rownames(newdat)[i]))
    if((length(num_1)==1)||(length(num_2)==1)){
        if(length(num_1)>0){
            nominator<-as.numeric(newintron[num_1,])
            denominator<-(as.numeric(newdat[i,])+as.numeric(newintron[num_1,]))
            tmp<-paste(nominator,denominator,sep=":")
        }
        if(length(num_2)>0){
            nominator<-as.numeric(newintron[num_2,])
            denominator<-(as.numeric(newdat[i,])+as.numeric(newintron[num_2,]))
            tmp<-paste(nominator,denominator,sep=":")
        }
        newdat[i,]<-tmp
    }
    if((length(num_1)==0)&&(length(num_2)==0)){
        newdat[i,1]<-NA
    }
}
newdat<-na.omit(newdat)

write.table(newdat,paste("/home/users/nus/e0950183/scratch/new_simulation_result/IR/sample400_iso2_cell",cell,"_effect_",effect,"/sample",sample,"_meta2_IR.site",sep=""),row.names=TRUE,col.names=TRUE,quote=FALSE)
