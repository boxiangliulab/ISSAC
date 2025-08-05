#!/bin/bash
#PBS -N Job_Name 
#PBS -q normal
#PBS -P 11003054
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=4:mem=32GB 
#PBS -o Out_File_Name_STAR.out 
#PBS -e Error_File_Name_STAR.err 

cd $PBS_O_WORKDIR

conda activate ~/ISSAC_env

effect="1_7"
sample_num=400
isoform=2
cell=100
mkdir /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/model/${effect}_${sample_num}_${cell}
mkdir /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/QTL/${effect}_${sample_num}_${cell}

/home/users/nus/e0950183/code/ISSAC_combine_tools_v2/QTL_mapping/build/ISSAC model \
       -s /home/users/nus/e0950183/scratch/new_simulation_result/metacell/pheno/sample${sample_num}_iso${isoform}_effect${effect}_meta2.filtered \
       -p /home/users/nus/e0950183/scratch/new_simulation_result/PC/effect${effect}_${sample_num}.PCs \
       -g /data/projects/11003054/e0950183/compare_site_intron_based/effect_size/PC/ISSAC_PC/GRM.txt \
       -n 400 \
       -u /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/model/${effect}_${sample_num}_${cell} \
       -o /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/model/model.middle


#list=/data/projects/11003054/e0950183/compare_QTLtools_GLMM/ISSAC_QTL/site.list

for i in `seq 1 22`;do
 ls /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/model/${effect}_${sample_num}_${cell}/*.middle | while read line;do 
  file_name=`ls $line | cut -d '/' -f 12`
  old_name=`echo $file_name | cut -d '.' -f 1`
  echo $old_name
  new_word="chr1:+:500"
  sed -i "s/$old_name/$new_word/g" "$line" #> /data/projects/11003054/e0950183/compare_QTLtools_GLMM/ISSAC_QTL/model_chr1/${file_name}
  mkdir /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/QTL/${effect}_${sample_num}_${cell}/${old_name}
  echo $old_name > test_${effect}_${sample_num}_${cell}.site
  /home/users/nus/e0950183/scratch/ISSAC_combine_tools_v2/QTL_mapping/build/ISSAC QTL -s test_${effect}_${sample_num}_${cell}.site \
       -p /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/model/${effect}_${sample_num}_${cell} \
       -m /home/users/nus/e0950183/scratch/new_simulation_result/sample_list/sample${sample_num}_meta.list \
       -c chr${i} \
       -x /home/users/nus/e0950183/scratch/new_simulation_result/PC/effect${effect}_${sample_num}.PCs \
       -v /data/projects/11003054/e0950183/compare_site_intron_based/effect_size/genotype/pseudobulk_geno.bcf \
       -w 1000 \
       -t 1 \
       -o /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/QTL/${effect}_${sample_num}_${cell}/${old_name}
 done
done



effect="1_7"
sample_num=400
isoform=2
cell=100
ISSAC=/home/users/nus/e0950183/code/ISSAC_combine_tools_v2/QTL_mapping/build/ISSAC
rm -rf /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/model/${effect}_${sample_num}_${cell}_IR
rm -rf /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/QTL/${effect}_${sample_num}_${cell}_IR

mkdir /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/model/${effect}_${sample_num}_${cell}_IR
mkdir /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/QTL/${effect}_${sample_num}_${cell}_IR

output_prefix=/home/users/nus/e0950183/scratch/new_simulation_result/IR/sample400_iso2_cell${cell}_effect_${effect}/sample${sample_num}_meta2_IR
out_file=${output_prefix}.filtered
prop_file=${output_prefix}.prop
sd=0.1
na_prop=0.4

splice_read_count=${output_prefix}.site
$ISSAC pheno_output $splice_read_count $out_file $prop_file $sd $na_prop

$ISSAC model \
       -s $out_file \
       -p /home/users/nus/e0950183/scratch/new_simulation_result/PC/effect${effect}_${sample_num}.PCs \
       -g /data/projects/11003054/e0950183/compare_site_intron_based/effect_size/PC/ISSAC_PC/GRM.txt \
       -n 400 \
       -u /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/model/${effect}_${sample_num}_${cell}_IR \
       -o /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/model/model.middle


#list=/data/projects/11003054/e0950183/compare_QTLtools_GLMM/ISSAC_QTL/site.list

for i in `seq 1 22`;do
 ls /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/model/${effect}_${sample_num}_${cell}_IR/*.middle | while read line;do 
  file_name=`ls $line | cut -d '/' -f 12`
  old_name=`echo $file_name | cut -d '.' -f 1`
  echo $old_name
  new_word="chr1:+:500"
  sed -i "s/$old_name/$new_word/g" "$line" #> /data/projects/11003054/e0950183/compare_QTLtools_GLMM/ISSAC_QTL/model_chr1/${file_name}
  mkdir /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/QTL/${effect}_${sample_num}_${cell}_IR/${old_name}
  echo $old_name > test_${effect}_${sample_num}_${cell}_IR.site
  /home/users/nus/e0950183/scratch/ISSAC_combine_tools_v2/QTL_mapping/build/ISSAC QTL -s test_${effect}_${sample_num}_${cell}_IR.site \
       -p /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/model/${effect}_${sample_num}_${cell}_IR \
       -m /home/users/nus/e0950183/scratch/new_simulation_result/sample_list/sample${sample_num}_meta.list \
       -c chr${i} \
       -x /home/users/nus/e0950183/scratch/new_simulation_result/PC/effect${effect}_${sample_num}.PCs \
       -v /data/projects/11003054/e0950183/compare_site_intron_based/effect_size/genotype/pseudobulk_geno.bcf \
       -w 1000 \
       -t 1 \
       -o /home/users/nus/e0950183/scratch/new_simulation_result/result/ISSAC/QTL/${effect}_${sample_num}_${cell}_IR/${old_name}
 done
done
