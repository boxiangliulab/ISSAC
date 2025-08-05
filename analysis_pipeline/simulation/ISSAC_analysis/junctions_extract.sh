#!/bin/bash
#PBS -N Job_Name 
#PBS -q normal
#PBS -P 11003054
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=2:mem=16GB 
#PBS -o Out_File_Name_STAR.out 
#PBS -e Error_File_Name_STAR.err 

cd $PBS_O_WORKDIR

samtools=~/ISSAC_env/bin/samtools
ISSAC=/home/users/nus/e0950183/code/ISSAC_combine_tools_v2/QTL_mapping/build/ISSAC


isoform=2
effect="1_7"
meta=2
cell_num=100
mkdir /home/users/nus/e0950183/scratch/new_simulation_result/metacell/iso${isoform}_meta${meta}/${effect}
mkdir /home/users/nus/e0950183/scratch/new_simulation_result/junc/iso${isoform}/${effect}_cell${cell_num}
#mkdir /data/projects/11003054/e0950183/compare_site_intron_based/effect_size/junc/metacell_junc/${effect}
for i in `seq 1 400`;do
   bamfile=/home/users/nus/e0950183/scratch/new_simulation_result/bam/iso2/${effect}_cell${cell_num}/S${i}Aligned.sortedByCoord.out.bam
   output_file=/home/users/nus/e0950183/scratch/new_simulation_result/junc/iso${isoform}/${effect}_cell${cell_num}/S${i}.junc
   ${samtools} index ${bamfile}
   ${ISSAC} junctools extract -a 8 -m 50 -s FR ${bamfile} -o ${output_file}
   barcode=/home/users/nus/e0950183/scratch/new_simulation_result/metacell/meta${meta}_1.list
   output=/home/users/nus/e0950183/scratch/new_simulation_result/metacell/iso${isoform}_meta${meta}/${effect}/S${i}:1.stat
   ${ISSAC} juncstat ${barcode} ${output_file} ${output}
   barcode=/home/users/nus/e0950183/scratch/new_simulation_result/metacell/meta${meta}_2.list
   output=/home/users/nus/e0950183/scratch/new_simulation_result/metacell/iso${isoform}_meta${meta}/${effect}/S${i}:2.stat
   ${ISSAC} juncstat ${barcode} ${output_file} ${output}
done
