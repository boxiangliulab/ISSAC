# ISSAC （Integrative single-cell splicing analysis and QTL caller） 

### ISSAC is designed for sQTL mapping at the single cell level by modeling metacell splice site usage ratio with generalized binomial mixed model.
![image](https://github.com/tibettiger/ISSAC/blob/ISSAC_analysis/img/figure1_v2.png)

## tutorial
## before use
Several C++ libraries needed to be added to your conda environment firstly before entering ISSAC world
1) htslib1.3
2) gsl
3) eigen3
4) nlopt
5) crypto

## How to install ISSAC?
## Metacell detection
To improve statistical power while maintaining low sparsity intrinsic to single-cell data, we prepared phenotype for downstream cis-sQTL based on the concept of metacells. For each donor, initial estimates for each cell cluster were obtained based on PC embedding obtained from snRNA-seq gene expression matrix through k nearest neighbors’ graph construction and Louvain clustering (with m as the default parameters of the size of each cluster). Then we iteratively decomposed clusters with less than m cells to adjacent larger clusters to ensure each cluster having at least m cells. Finally, each cluster with at least m cells was treated as one metacell for downstream sQTL mapping tasks. The original code of our metacell detection method is uploaded: https://github.com/tibettiger/ISSAC/blob/ISSAC_analysis/analysis_pipeline/DLPFC_analysis/sQTL_mapping/metacell_calling.py 

## Phenotype preparation
After obtaining bam files for each metacell, we first perform junctions extract and quantification.

1）
```
function test() {
  console.log("notice the blank line before this function?");
}
```

`junctools=junctools`
juncstat=junc_stat  
bamfile=*.bam  //metacell level bamfile; needing .bai index file
output_file=*.junc
${junctools} junctions extract -a 8 -m 50 -s FR ${bamfile} -o ${output_file} (FR for 5' scRNA-seq & RF for 3' scRNA-seq)
barcode=*.barcode //the cell barcode you want to extract
output=*.stat
${juncstat} ${barcode} ${output_file} ${output}`

Examples of junctools output
Each column‘s meaning: 1) chrom 2)start of the read 3) end of the read 4) junc seq 5) read count
6) strand 7) start of the junction 8) end of the junction 9,10) no usage 11) anchor length of the read
12) no usage 13) cell barcode of the read 14) UMI information of the read

Examples of stat output
Each column‘s meaning: 1) intron junction information 2) UMI count of the intron in the bam files
Here, we use UMI counts to quantify intron junction reads to avoid PCR amplication bias

After obtaining all the .junc and .stat file for each metacell, we perform intron cluster and site quantification

2)
phenotype_group=phenotype_group
read_threshold=30 (could be adjust; intron with its whole supported UMI counts less than the threshold will be deleted)
junc_pos=/${}_junc (the position of all the .stat file)
sample_file=${}_sample (deposit the sample name)
out_file_prefix=${} (the output file’s prefix)
log_out=${} (the log file)
$phenotype_group $sample_file $junc_pos $out_file_prefix $read_threshold $log_out

Output of 03_phenotype_group include 4 files;
.inclu_exclu (deposit each splice site’s info including included intron and excluded intron)
.intron.out (deposit intron cluster’s results and UMI counts in each sample)
.refined (deposit each intron’s total number larger than the threshold)
.site (deposit the phenotype of each splice site)

After obtaining initial phenotype for each site (.site), we perform phenotype filtering 
3)
splice_read_count=*.site  (the initial phenotype file obtained from previous step)
phenotype_output=phenotype_output

out_file=*.filtered (deposit filtered phenotype file)
prop_file=*.prop
sd=0.01 (splice sites with variations less than the threshold will be deleted)
na_prop=0.4 (splice sites lacking phenotypes in larger than the threshold of the whole samples will be deleted)

$phenotype_output $splice_read_count $out_file $prop_file $sd $na_prop

The output of filter including two files:
(1).filtered (2).prop

## model construction
After obtaining phenotype(.filtered), GRM, genotype(.bcf &.csi),  .PC file, we could perform model construction for each splice site to estimate coefficients of fixed effect and random effect
GRM preparation

use plink to obtain .grm.bin, .grm.grm.N.bin file (plink --bfile pruned_output --make-grm-bin --out grm_output)
use below file to obtain GRM matrix 

Examples of GRM file: Note: all the phenotypes and PCs’ names should be composed of ${sample}:${seq}; and ${sample} must exist in GRM file

ISSAC=ISSAC
$ISSAC model –s *.filtered   phenotype    
       -p *.PC    PC file  
       -g GRM.txt    GRM file  
       -n 617   sample numbers of GRM file  
       -u ${model_pos}  the position of all the model files   
       -o  middle  the output name  


Note: ISSAC model will output all the splice site’s model containing in the phenotype file

Example of model files output by ISSAC model:

The first line includes: 1) splice site 2) dispersion parameter 3) variance components of random effects
another four lines are 1) residuals 2) pi 3) total UMI count for the site in each sample 4) used UMI count for the site in each sample

## score tests for sQTL mapping

After obtaining model files, then we could utilize bcf file to perform cis-sQTL mapping

ISSAC QTL mapping:

$ISSAC QTL   
       -p $model_file_pos   the model files position  
       -m  *.common the sample names of phenotype and PC  file  
       -c chr${i}  the chromosome of splice sites  
       -x *.PC   PC file  
       -v chr${i}.recode.bcf  genotype file  
       -w 1000000  the windows within this range will be used to perform cis-sQTL mapping  
       -o $output_pos   the output file’s position  
       -t 1  the sQTLs with pvalue less than this threshold will be output to result files  

Output of ISSAC QTL includes .result file and it will output results of all the splice sites in your desired chromosome of your .site file

Example of result file: 1) splice site 2) SNP 3) pvalue 4) effect size 5) standard error of effect size

Then the associations between the desired splice sites and the corresponding SNPs are determined.










