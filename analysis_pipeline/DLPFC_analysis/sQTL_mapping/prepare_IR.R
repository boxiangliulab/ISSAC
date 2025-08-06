cell<-"astrocytes"

library(data.table)

detect_site<-c()
splice_site<-read.table(paste("/home/users/nus/e0950183/scratch/brain_major_sub/phenotype/",cell,"_single_intron_site",sep=""))
sample_name<-read.table(paste("/home/users/nus/e0950183/scratch/brain_major_sub/phenotype/",cell,"_sample",sep=""))

newdat<-as.data.frame(matrix(0,nrow(splice_site),nrow(sample_name)))
rownames(newdat)<-splice_site$V1

#intron<-fread(paste("/home/users/nus/e0950183/scratch/brain_major_sub/phenotype/",cell,".intron.out",sep=""))


#colnames(intron)[3:(nrow(sample_name)+2)]<-sample_name$V1

colnames(newdat)<-sample_name$V1

barcode<-read.table(paste("/data/projects/11003054/e0950183/brain_sQTL/metacell/",cell,"_meta.txt",sep=""),header=TRUE)
barcode$batch<-lapply(barcode$batch_ind, function(x) strsplit(x, ":")[[1]][1])
barcode$ind<-lapply(barcode$batch_ind, function(x) strsplit(x, ":")[[1]][2])

barcode$lib<-NA
barcode$tmp_barcode<-NA
barcode$new_barcode<-NA

num<-which(barcode$batch=="Source")
barcode$lib[num]<-lapply(barcode$barcode[num], function(x) strsplit(x, "_")[[1]][1])
barcode$tmp_barcode[num]<-as.character(lapply(barcode$barcode[num], function(x) strsplit(x, "_")[[1]][2]))
barcode$new_barcode[num]<-as.character(lapply(barcode$tmp_barcode[num], function(x) strsplit(x, "-")[[1]][1]))
barcode$tmp_barcode[num]<-paste(barcode$ind[num],barcode$batch[num],barcode$lib[num],sep="_")

num<-which(barcode$batch=="Target")
barcode$new_barcode[num]<-as.character(lapply(barcode$barcode[num], function(x) strsplit(x, "-")[[1]][1]))

##add source
num<-which(barcode$batch=="Source")
sample_name$ind<-lapply(sample_name$V1, function(x) strsplit(x, ":")[[1]][1])

###NG
nonsplit_sample<-unique(barcode$tmp_barcode[num])
###read in nonsplit
for(i in 1:length(nonsplit_sample)){
    print(i)
    nonsplit<-fread(paste("/home/users/nus/e0950183/scratch/brain_major_sub/phenotype/non_split_NG/",cell,"/",nonsplit_sample[i],".nonsplit",sep=""),sep=":")
    nonsplit<-as.data.frame(nonsplit)
    colnames(nonsplit)<-c("chr","pos","strand","barcode","other")
    nonsplit$site<-paste(nonsplit$chr,nonsplit$strand,nonsplit$pos,sep=":")
    ind<-strsplit(nonsplit_sample[i],"_")[[1]][1]
    lib<-strsplit(nonsplit_sample[i],"_")[[1]][3]
    num_barcode<-intersect(which(barcode$lib==lib),which(barcode$ind==ind))
    metacell<-unique(barcode$identify[num_barcode])
    for(j in 1:length(metacell)){
        num_metacell<-intersect(which(barcode$identify==metacell[j]),which(barcode$lib==lib))
        barcode_meta<-barcode$new_barcode[num_metacell]
        num_metacell<-which(nonsplit$barcode%in%barcode_meta)
        if(length(num_barcode)>0){
        new_nonsplit<-nonsplit[num_metacell,]
        site_metacell<-unique(nonsplit$site[num_metacell])
        detect_site<-unique(c(detect_site,site_metacell))
        counts<-as.data.frame(table(new_nonsplit$site))
        newdat[as.character(counts$Var1),metacell[j]]<-newdat[as.character(counts$Var1),metacell[j]]+counts$Freq
        }
    }
}
##add Target
num<-which(barcode$batch=="Target")
ind<-as.character(unique(barcode$ind[num]))
file<-list.files(paste("/home/users/nus/e0950183/scratch/brain_major_sub/phenotype/non_split_Cell/",cell,sep=""))
newfile<-as.data.frame(matrix(NA,length(file),2))
newfile$V1<-file
newfile$V2<-as.character(lapply(newfile$V1, function(x) strsplit(x, "_")[[1]][1]))

for(i in 1:length(ind)){
    print(ind[i])
    num<-which(newfile$V2==ind[i])
    if(length(num)>0){
    tmp<-as.data.frame(matrix(NA,0,5))
    for(j in 1:length(num)){
        print(j)
    nonsplit<-fread(paste("/home/users/nus/e0950183/scratch/brain_major_sub/phenotype/non_split_Cell/",cell,"/",newfile$V1[num[j]],sep=""),sep=":")
    nonsplit<-as.data.frame(nonsplit)
    tmp<-rbind(tmp,nonsplit)
    }
    nonsplit<-tmp
    #nonsplit$barcode<-as.character(lapply(nonsplit$V1, function(x) strsplit(x, ":")[[1]][4]))
    #nonsplit$chr<-as.character(lapply(nonsplit$V1, function(x) strsplit(x, ":")[[1]][1]))
    #nonsplit$strand<-as.character(lapply(nonsplit$V1, function(x) strsplit(x, ":")[[1]][3]))
    #nonsplit$pos<-as.character(lapply(nonsplit$V1, function(x) strsplit(x, ":")[[1]][2]))
    colnames(nonsplit)<-c("chr","pos","strand","barcode","other")
    nonsplit$site<-paste(nonsplit$chr,nonsplit$strand,nonsplit$pos,sep=":")
    num_barcode<-intersect(which(barcode$batch=="Target"),which(barcode$ind==ind[i]))
    metacell<-unique(barcode$identify[num_barcode])
    for(k in 1:length(metacell)){
        num_metacell<-which(barcode$identify==metacell[k])
        barcode_meta<-barcode$new_barcode[num_metacell]
        num_metacell<-which(nonsplit$barcode%in%barcode_meta)
        new_nonsplit<-nonsplit[num_metacell,]
        site_metacell<-unique(nonsplit$site[num_metacell])
        detect_site<-unique(c(detect_site,site_metacell))
        if(length(site_metacell)>0){
            print("contain")
            counts<-as.data.frame(table(new_nonsplit$site))
            newdat[as.character(counts$Var1),metacell[k]]<-newdat[as.character(counts$Var1),metacell[k]]+counts$Freq
        }
    }

    }
}



newdat<-newdat[detect_site,]


write.table(newdat,paste("/home/users/nus/e0950183/scratch/brain_major_sub/phenotype/",cell,".nonsplit",sep=""),row.names=TRUE,col.names=TRUE,quote=FALSE)
