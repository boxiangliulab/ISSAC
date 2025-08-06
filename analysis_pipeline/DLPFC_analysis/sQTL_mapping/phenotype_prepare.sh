#!/bin/bash
#PBS -N Job_Name 
#PBS -q large
#PBS -P 11003054
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=12:mem=1028GB 
#PBS -o Out_File_Name.out 
#PBS -e Error_File_Name.err 

cd $PBS_O_WORKDIR

cell="vascular"
##phenotype_group=/data/projects/11003054/e0950183/code/03_phenotype_group
ISSAC=/home/users/nus/e0950183/code/ISSAC_combine_tools_v2/QTL_mapping/build/ISSAC

ls /home/users/nus/e0950183/scratch/brain_major_sub/stat/major/celltype/${cell}/*.stat | cut -d '/' -f 12 | cut -d '.' -f 1 > \
  /home/users/nus/e0950183/scratch/brain_major_sub/phenotype/${cell}_sample

output_prefix=/home/users/nus/e0950183/scratch/brain_major_sub/phenotype/${cell}

log_out=${cell}.log

read_threshold=20

junc_pos=/home/users/nus/e0950183/scratch/brain_major_sub/stat/major/celltype/${cell}
sample_file=/home/users/nus/e0950183/scratch/brain_major_sub/phenotype/${cell}_sample
$ISSAC pheno_group $sample_file $junc_pos $output_prefix $read_threshold $log_out

out_file=${output_prefix}.filtered
prop_file=${output_prefix}.prop
sd=0.001
na_prop=0.5

splice_read_count=${output_prefix}.site

$ISSAC pheno_output $splice_read_count $out_file $prop_file $sd $na_prop



####IR
ISSAC=/home/users/nus/e0950183/code/ISSAC_combine_tools_v2/QTL_mapping/build/ISSAC
cell="oligodendrocyte"
output_prefix=/home/users/nus/e0950183/scratch/brain_major_sub/phenotype/${cell}_IR
out_file=${output_prefix}.filtered
prop_file=${output_prefix}.prop
sd=0.001
na_prop=0.5

splice_read_count=${output_prefix}.site

$ISSAC pheno_output $splice_read_count $out_file $prop_file $sd $na_prop
