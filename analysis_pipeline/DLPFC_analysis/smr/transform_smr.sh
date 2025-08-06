###perform plink to vcf

plink=/home/e0950183/miniconda3/envs/py27/bin/plink

ls /home/e0950183/project/brain_sQTL/coloc/smr/vcf/*.vcf | while read line;do
   file_name=`echo $line | cut -d '.' -f 1`
   $plink --vcf $line --make-bed --out $file_name
done


###

smr=/home/e0950183/project/brain_sQTL/coloc/smr/smr_linux_x86_64 

ls /home/e0950183/project/brain_sQTL/coloc/smr/smr_sQTL/* | while read line;do
  $smr --qfile $line \
   --make-besd --out $line
done 


cat /home/e0950183/project/brain_sQTL/coloc/smr/smr_list | while read line;do
 chr1=`echo $line | cut -d ':' -f 1`
 chr=`echo $chr1 | cut -d '-' -f 2`
 pos1=`echo $line | cut -d ':' -f 3`
 pos=`echo $pos1 | cut -d '-' -f 1` 
 $smr --bfile /home/e0950183/project/brain_sQTL/coloc/smr/vcf/${chr}_${pos} \
   --gwas-summary /home/e0950183/project/brain_sQTL/coloc/smr/smr_GWAS/$line \
   --beqtl-summary /home/e0950183/project/brain_sQTL/coloc/smr/smr_sQTL/$line \
   --out /home/e0950183/project/brain_sQTL/coloc/smr/result/${line} --thread-num 10 --diff-freq-prop 1 --peqtl-smr 1e-3
done


cat /home/e0950183/project/brain_sQTL/coloc/summary_coloc_H4_over_0_75.txt | awk '{print $14 "-" $2 "-" $13}' | while read line;do
   num=`ls /home/e0950183/project/brain_sQTL/coloc/smr/smr_GWAS/* | grep $line | wc -l`
   if [ $num -eq 0 ];then
     echo $line
   fi
done
