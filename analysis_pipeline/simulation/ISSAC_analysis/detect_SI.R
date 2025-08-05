###Shell code to extract nonsplit reads crossing one splice site
#!/bin/bash
#PBS -N Job_Name 
#PBS -q normal
#PBS -P 11003054
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=2:mem=16GB 
#PBS -o Out_File_Name_STAR.out 
#PBS -e Error_File_Name_STAR.err 

cd $PBS_O_WORKDIR


ISSAC=/home/users/nus/e0950183/code/ISSAC_combine_tools_v2/QTL_mapping/build/ISSAC
###

iso=2
sample=400
cell=100
effect="1_5"
site=/home/users/nus/e0950183/scratch/new_simulation_result/IR/sample${sample}_iso${iso}_cell${cell}_effect_${effect}/single_intron_site
mkdir /home/users/nus/e0950183/scratch/new_simulation_result/IR/sample${sample}_iso${iso}_cell${cell}_effect_${effect}
barcode=/home/users/nus/e0950183/scratch/new_simulation_result/metacell/synthetic_cell_barcode.txt

for i in `seq 1 400`;do
  line=/data/projects/11003054/e0950183/compare_site_intron_based/effect_size/bam/${effect}/S${i}Aligned.sortedByCoord.out.bam
  #line=/home/users/nus/e0950183/scratch/new_simulation_result/bam/iso2/${effect}_cell${cell}/S${i}Aligned.sortedByCoord.out.bam
  output_file=/home/users/nus/e0950183/scratch/new_simulation_result/IR/sample${sample}_iso${iso}_cell${cell}_effect_${effect}/S${i}.nonsplit
  $ISSAC IR extract -s FR -b ${barcode} -t ${site} ${line} -o ${output_file}
done






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
