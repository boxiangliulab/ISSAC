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




