cell="microglia"

mkdir /home/users/nus/e0950183/scratch/brain_major_sub/phenotype/non_split_NG/${cell}

for i in `seq 0 42`;do
   start=$((i*10+1))
   end=$((i*10+10))
   echo -e "
#!/bin/bash
#PBS -N Job_Name 
#PBS -q normal
#PBS -P 11003054
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=4:mem=32GB 
#PBS -o Out_File_Name_STAR.out 
#PBS -e Error_File_Name_STAR.err 

cd \$PBS_O_WORKDIR

ISSAC=/home/users/nus/e0950183/code/ISSAC_combine_tools_v2/QTL_mapping/build/ISSAC
samtools=~/ISSAC_env/bin/samtools
site=/home/users/nus/e0950183/scratch/brain_major_sub/phenotype/${cell}_single_intron_site

cat /data/projects/11003054/e0950183/brain_sQTL/metacell/${cell}_meta.txt | grep "Source" | awk '{print \$2}' | cut -d ':' -f 2 | sort | uniq | sed -n '${start},${end}p' | while read sample;do
  cat /data/projects/11003054/e0950183/brain_sQTL/metacell/${cell}_meta.txt | grep "Source" | grep \$sample | awk -F '-1' '{print \$1}' > \
     /home/users/nus/e0950183/scratch/brain_major_sub/phenotype/non_split/barcode/${cell}_\${sample}_Source.barcode
  ls /home/users/nus/e0950183/scratch/NG_bam/*\${sample}* | while read bam;do
     \$samtools index \$bam
     bam_prefix=\`echo \$bam | cut -d '/' -f 8 | cut -d '_' -f 1\`
     echo \$bam_prefix
     cat /home/users/nus/e0950183/scratch/brain_major_sub/phenotype/non_split/barcode/${cell}_\${sample}_Source.barcode | grep \$bam_prefix | cut -d '_' -f 2 > \
       /home/users/nus/e0950183/scratch/brain_major_sub/phenotype/non_split/barcode/${cell}_\${sample}_\${bam_prefix}_Source.barcode
     barcode=/home/users/nus/e0950183/scratch/brain_major_sub/phenotype/non_split/barcode/${cell}_\${sample}_\${bam_prefix}_Source.barcode
     \$ISSAC IR extract -s RF -b \${barcode} -t \${site} \${bam} -o \
        /home/users/nus/e0950183/scratch/brain_major_sub/phenotype/non_split_NG/${cell}/\${sample}_Source_\${bam_prefix}.nonsplit
  done
done" > /home/users/nus/e0950183/scratch/batch_scripts/junc_stat_Source_${cell}_${i}.sh
   qsub /home/users/nus/e0950183/scratch/batch_scripts/junc_stat_Source_${cell}_${i}.sh
done
