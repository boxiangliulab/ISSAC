#!/bin/bash
celltype="OPC"
rm -rf /home/users/nus/e0950183/scratch/brain_major_sub/result/major/${celltype}_model
mkdir /home/users/nus/e0950183/scratch/brain_major_sub/result/major/${celltype}_model
for i in `seq 1 1`;do
  echo -e "
#!/bin/bash
#PBS -N Job_Name 
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=4:mem=32GB 
#PBS -P 11003054
#PBS -o Out_File_Name.out 
#PBS -e Error_File_Name.err 

cd $PBS_O_WORKDIR
conda activate /home/users/nus/e0950183/ISSAC_env
/home/users/nus/e0950183/code/ISSAC_combine_tools_v2/QTL_mapping/build/ISSAC model \
       -s /home/users/nus/e0950183/scratch/brain_major_sub/phenotype/${celltype}.filtered \
       -p /data/projects/11003054/e0950183/brain_sQTL/PC/PC/${celltype}.PC \
       -g /data/projects/11003054/e0950183/brain_sQTL/GRM/GRM_579.txt \
       -n 579 \
       -u /home/users/nus/e0950183/scratch/brain_major_sub/result/major/${celltype}_model \
       -o /home/users/nus/e0950183/scratch/brain_major_sub/result/major/${celltype}_model
" > /home/users/nus/e0950183/scratch/batch_scripts/run_split_model_${celltype}_${i}.sh
   qsub /home/users/nus/e0950183/scratch/batch_scripts/run_split_model_${celltype}_${i}.sh
done
