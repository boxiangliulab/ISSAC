#!/bin/bash
#PBS -N Job_Name 
#PBS -q large
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=108GB 
#PBS -o Out_File_Name_STAR.out 
#PBS -e Error_File_Name_STAR.err 

cd $PBS_O_WORKDIR

STAR=~/ISSAC_env/bin/STAR
genomeDir=/data/projects/11003054/share/data/reference/STAR_index/hg38_gencode_v32
gtfFile=/data/projects/11003054/share/data/reference/hg38/gencode.v32.primary_assembly.annotation.gtf  
whitelist=/data/projects/11003054/e0950183/compare_site_intron_based/JP_RIK_H002_CD4+_T_naive_simulation/tmp_scisosim/synthetic_cell_barcode.txt

iso=2
effect="1_1"
mkdir /data/projects/11003054/e0950183/compare_site_intron_based/effect_size/bam/${effect}
for i in `seq 301 400`;do
   readFilesIn1=/data/projects/11003054/e0950183/compare_site_intron_based/effect_size/fq/${effect}/Isoform${i}.read1.bed2fa.sorted.fq
   readFilesIn2=/data/projects/11003054/e0950183/compare_site_intron_based/effect_size/fq/${effect}/Isoform${i}.read2.bed2fa.sorted.fq
   $STAR \
    --runThreadN 8 \
    --genomeDir  $genomeDir \
    --twopassMode Basic \
    --soloStrand Forward \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist $whitelist --soloCBmatchWLtype 1MM --soloUMIdedup 1MM_Directional_UMItools \
    --outSAMstrandField intronMotif \
    --soloBarcodeMate 1 --clip5pNbases 39 0 \
    --readFilesIn $readFilesIn1 $readFilesIn2  \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ XS\
    --sjdbGTFfile $gtfFile \
    --soloCBstart 1 \
    --soloCBlen 16 \
    --soloUMIstart 17 \
    --soloUMIlen 10 \
    --outFileNamePrefix /data/projects/11003054/e0950183/compare_site_intron_based/effect_size/bam/${effect}/S${i}  
done
