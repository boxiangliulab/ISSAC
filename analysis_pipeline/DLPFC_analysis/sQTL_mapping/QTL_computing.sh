celltype="microglia"

mkdir -p batch_scripts
rm -rf /data/projects/11003054/e0950183/brain_sQTL/QTL_result/${celltype}
mkdir /home/users/nus/e0950183/scratch/brain_major_sub/result/major/${celltype}_qtl
ls /home/users/nus/e0950183/scratch/brain_major_sub/result/major/${celltype}_model/*.middle | \
  cut -d '.' -f 1 | cut -d '/' -f 11 > /home/users/nus/e0950183/scratch/brain_major_sub/result/major/${celltype}.site

cat /home/users/nus/e0950183/scratch/brain_major_sub/phenotype/${celltype}.newfiltered | head -n 1 | tr ' ' '\n' > \
   /home/users/nus/e0950183/scratch/brain_major_sub/phenotype/${celltype}.common
for i in `seq 1 22`;do
    line=/home/users/nus/e0950183/scratch/brain_major_sub/result/major/${celltype}_lack_site
    echo -e "
#!/bin/bash
#PBS -N Job_Name 
#PBS -q normal
#PBS -P 11003054
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=4:mem=32GB 
#PBS -o Out_File_Name.out 
#PBS -e Error_File_Name.err 

cd \$PBS_O_WORKDIR
conda activate /home/users/nus/e0950183/ISSAC_env
/home/users/nus/e0950183/code/ISSAC/build/ISSAC QTL -s ${line} \
       -p /home/users/nus/e0950183/scratch/brain_major_sub/result/major/${celltype}_model \
       -m /home/users/nus/e0950183/scratch/brain_major_sub/phenotype/${celltype}.common \
       -x /home/users/nus/e0950183/scratch/brain_major_sub/PC/${celltype}.PC \
       -c chr${i} \
       -v /data/projects/11003054/e0950183/brain_sQTL/genotype/after_imputation_qc/chr${i}.filtered.bcf \
       -w 1000000 \
       -o /home/users/nus/e0950183/scratch/brain_major_sub/result/major/${celltype}_qtl \
       -t 1 
" > /home/users/nus/e0950183/scratch/batch_scripts/QTL_${celltype}_chr${i}.sh
    qsub /home/users/nus/e0950183/scratch/batch_scripts/QTL_${celltype}_chr${i}.sh
done
