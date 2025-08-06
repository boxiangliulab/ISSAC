cell="astrocytes"
mkdir /home/users/nus/e0950183/scratch/brain_major_sub/stat/major/barcode/${cell}
mkdir /home/users/nus/e0950183/scratch/brain_major_sub/stat/major/celltype/${cell}

for i in `seq 0 18`;do
   start=$((i*100+1))
   end=$((i*100+100))
   echo -e "
#!/bin/bash
#PBS -N Job_Name 
#PBS -P 11003054
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=2:mem=16GB 
#PBS -o Out_File_Name.out 
#PBS -e Error_File_Name.err 

cd \$PBS_O_WORKDIR

conda activate ~/ISSAC_env
ISSAC=/data/projects/11003054/e0950183/code/ISSAC

cat /data/projects/11003054/e0950183/brain_sQTL/metacell/${cell}_meta.txt | grep "Source" | awk '{print \$4}' | sort | uniq | sed -n '${start},${end}p' | while read line;do
  cat /data/projects/11003054/e0950183/brain_sQTL/metacell/${cell}_meta.txt | awk -v line="\$line" '{if(\$4==line)print \$1}' | sed 's/...$//' > \
    /home/users/nus/e0950183/scratch/brain_major_sub/stat/major/barcode/${cell}/\${line}.barcode
  barcode=/home/users/nus/e0950183/scratch/brain_major_sub/stat/major/barcode/${cell}/\${line}.barcode
  sample=\`echo \$line | cut -d ':' -f 1\`
  output_file=/data/projects/11003054/e0950183/brain_sQTL/junc/total/NG/\${sample}_batch.junc
  output=/home/users/nus/e0950183/scratch/brain_major_sub/stat/major/celltype/${cell}/\${line}.stat
  \${ISSAC} juncstat \${barcode} \${output_file} \${output}
done" > /home/users/nus/e0950183/scratch/batch_scripts/junc_stat_${cell}_${i}.sh
  qsub /home/users/nus/e0950183/scratch/batch_scripts/junc_stat_${cell}_${i}.sh
done
