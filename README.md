# ISSAC （Integrative single-cell splicing analysis and QTL caller） 

### ISSAC is a scalable tool for single-cell sQTL mapping by modeling metacell splice site usage ratio with generalized binomial mixed model.
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
After obtaining bam files for each metacell, we first perform junctions extract.

```
bamfile=*.bam  //bamfile; needing .bai index file
output_file=*.junc
${ISSAC} junctions extract -a 8 -m 50 -s FR ${bamfile} -o ${output_file} (FR for 5' scRNA-seq & RF for 3' scRNA-seq)
barcode=*.barcode //the cell barcode you want to extract
output=*.stat
${ISSAC} juncstat ${barcode} ${output_file} ${output}`
```

Examples of junctions output:
```
chr1    6186734 6193054 JUNC00000065    4       -       6186816 6192929 255,0,0 2       82,125  0,6195  TATTACCCAGCCAATT        TGCTAGCGCT
chr1    6186737 6193054 JUNC00000066    5       -       6186816 6192929 255,0,0 2       79,125  0,6192  ACCCACTGTAATAGCA        ACTGAAATTT
chr1    6186743 6193007 JUNC00000069    2       -       6186816 6192929 255,0,0 2       73,78   0,6186  TATTACCCAGCCAATT        GTCATGTAAT
chr1    6186743 6193041 JUNC00000067    3       -       6186816 6192929 255,0,0 2       73,112  0,6186  CGGGTCACAACGATGG        TTAGTTCCCT
chr1    6186743 6193054 JUNC00000068    4       -       6186816 6192929 255,0,0 2       73,125  0,6186  CGGGTCACAACGATGG        TCGAGAGTAC
chr1    6186745 6193009 JUNC00000070    1       -       6186816 6192929 255,0,0 2       71,80   0,6184  TATTACCCAGCCAATT        AAGGTACAAC
chr1    6186751 6193019 JUNC00000071    3       -       6186816 6192929 255,0,0 2       65,90   0,6178  ACCCACTGTAATAGCA        GGGATTCAAT
chr1    6186751 6193019 JUNC00000073    4       -       6186816 6192929 255,0,0 2       65,90   0,6178  CGGGTCACAACGATGG        GCGCAACAAT
chr1    6186751 6193036 JUNC00000072    2       -       6186816 6192929 255,0,0 2       65,107  0,6178  ACCCACTGTAATAGCA        GGGAATCAGT
chr1    6186754 6193018 JUNC00000075    1       -       6186816 6192929 255,0,0 2       62,89   0,6175  ACCCACTGTAATAGCA        ATCACCTTTT
chr1    6186754 6193018 JUNC00000077    1       -       6186816 6192929 255,0,0 2       62,89   0,6175  TATTACCCAGCCAATT        TGCTATAACA
chr1    6186754 6193046 JUNC00000074    5       -       6186816 6192929 255,0,0 2       62,117  0,6175  ACCCACTGTAATAGCA        TCATTTATAT
chr1    6186754 6193046 JUNC00000076    3       -       6186816 6192929 255,0,0 2       62,117  0,6175  GCTGGGTTCGTTACAG        TGACCTACCG
chr1    6186754 6193054 JUNC00000078    4       -       6186816 6192929 255,0,0 2       62,125  0,6175  CGGGTCACAACGATGG        TCTTTCACAA
```
Each column‘s meaning: 1) chrom 2)start of the read 3) end of the read 4) junctions seq 5) numbers of junction read counts with the same UMI and cell barcode mapped to the same junction
6) strand 7) start of the junction 8) end of the junction 9,10) no usage 11) anchor length of the read
12) no usage 13) cell barcode of the read 14) UMI information of the read

After obtaining initial junctions information, we used juncstat to count CB-UMI based junction reads.
Examples of stat output:
```
chr10:+:124484326:124488421     1
chr10:+:124801901:124806830     1
chr10:+:124806921:124816575     1
chr10:+:124816612:124819383     1
chr10:+:12553356:12666735       1
chr10:+:125719915:125721203     3
chr10:+:130136512:130145210     1
chr10:+:131901155:131901237     1
chr10:+:131901330:131911561     2
chr10:+:132397351:132404625     2
chr10:+:132538171:132607914     2
chr10:+:132749612:132749770     1
chr10:+:132749845:132777670     1
chr10:+:132777782:132780848     1
```
Each column‘s meaning: 1) intron junction information 2) CB-UMI based counts of the intron in the bam files
Here, we use CB-UMI based counts to quantify intron junction reads to avoid PCR amplication bias

After obtaining all the .junc and .stat file for each metacell, we perform splice site partner definition and splice usage ratio computing

```
read_threshold=20 ## could be adjust; intron with its whole supported UMI counts less than the threshold will not be considered in the splice site partner definition
junc_pos=/directory/junc ## the directory of all the .stat file
sample_file=*_sample ## deposit the sample name
out_file_prefix=* ## the output file’s prefix
log_out=* # the log file of intermediate results
$ISSAC pheno_group $sample_file $junc_pos $out_file_prefix $read_threshold $log_out
```

Output of ISSAC pheno_group include 4 files;
.inclu_exclu (deposit each splice site’s info including included intron and excluded intron)  
The file of .inclu_exclu is shown as below:
```
chr10:+:73 110985627
included 110919657:110985627 
excluded 110919657:110951607 110919657:110964124 
chr10:+:76 111114416
included 111114416:111114902 111114416:111117959 
excluded 
chr10:+:76 111114902
included 111114416:111114902 
excluded 111114416:111117959 
chr10:+:76 111115030
included 111115030:111117959 
excluded 111114416:111117959 
chr10:+:76 111117959
included 111114416:111117959 111115030:111117959 
excluded 
chr10:+:79 11165682
included 11165682:11210670 11165682:11217424 11165682:11249152 
excluded 
chr10:+:79 11210670
included 11165682:11210670 
excluded 11165682:11217424 11165682:11249152 
chr10:+:79 11217424
included 11165682:11217424 
excluded 11165682:11249152 
chr10:+:79 11217507
included 11217507:11249152 
excluded 11165682:11249152 

```
Splice sites(first line), junction reads supporting the usage of the site (second line) & junction reads competing the usage of the site (third line) are involved in the file.
The file of .intron.out is shown as below:
```
chr10:+:231 chr10:+:22609717:22714182 0 0 0 0 0 0
chr10:+:232 chr10:+:22853655:22868088 0 0 0 0 0 0
chr10:+:233 chr10:+:22925636:22927862 0 0 0 0 0 0
chr10:+:233 chr10:+:22925636:22931995 2 0 2 0 1 0
chr10:+:233 chr10:+:22928106:22931995 0 0 0 0 0 0
chr10:+:234 chr10:+:22932044:22946143 1 1 1 0 0 0
chr10:+:234 chr10:+:22932044:22955806 0 0 0 0 0 0
chr10:+:234 chr10:+:22946261:22955806 1 0 0 0 0 0
chr10:+:235 chr10:+:22955932:22959069 0 0 0 0 0 0
chr10:+:236 chr10:+:22959138:22959398 0 0 0 0 0 0
chr10:+:237 chr10:+:23095726:23104143 0 0 0 0 0 0
chr10:+:238 chr10:+:23104244:23110241 0 0 1 0 0 0
```
.intron.out deposits intron clusters' results and UMI counts for each intron in each sample.
The file of .refined is shown as below:
```
chr10:+:111 11852215:11862911:125 11852215:11866530:640
chr10:+:112 118730410:118754395:128
chr10:+:113 119104197:119104763:91 119104197:119107947:44
chr10:+:114 119108164:119111848:35
chr10:+:115 119207969:119326515:530 119207969:119380814:38 119326611:119380814:182
chr10:+:116 119380927:119396694:105
chr10:+:117 119396772:119423165:72
chr10:+:118 119423266:119424992:58
chr10:+:119 119425085:119430374:58
chr10:+:12 101587457:101594657:27 101588400:101594657:263
```
.refined deposits each intron's CB-UMI based junction reads across all samples
The file of .site is shown as below:
```
chr10:-:14553321 0:0 1:1 0:1 0:0 0:0 0:0 0:0
chr10:-:14553387 0:0 1:1 0:1 0:0 0:0 0:0 0:0
chr10:-:16752538 0:0 0:0 0:0 0:0 0:0 0:0 0:0
chr10:-:16764510 0:0 0:0 1:1 0:0 0:0 0:0 0:0
chr10:-:16782084 13:13 11:11 17:18 9:9 6:6 15:15 15:15
chr10:-:16816972 11:13 11:11 14:18 8:9 6:6 15:15 14:15
chr10:-:16817084 9:11 12:12 13:17 8:9 7:7 13:13 15:16
chr10:-:17359600 0:0 0:0 0:0 0:0 0:0 0:0 0:0
chr10:-:17387100 0:0 0:0 0:0 0:0 0:0 0:0 0:0
chr10:-:19684254 0:1 0:0 0:0 0:0 0:3 0:2 0:3
chr10:-:19706080 0:0 0:0 0:0 0:0 0:0 0:0 0:0
chr10:-:19711288 0:0 0:0 0:0 0:0 0:0 0:0 1:1
chr10:-:19714521 0:0 0:0 0:0 0:0 1:1 0:0 0:1
chr10:-:19716890 1:1 0:0 0:0 0:0 2:3 2:2 2:3
chr10:-:19718499 0:1 0:0 0:0 0:0 1:3 2:2 2:3
chr10:-:19718506 1:1 0:0 0:0 0:0 2:3 0:2 1:3
chr10:-:19686434 0:0 0:0 0:0 0:0 0:0 0:0 0:0
chr10:-:19701526 0:0 0:0 0:0 0:0 0:0 0:0 0:0
chr10:-:19711288 1:1 0:0 0:0 0:0 0:2 2:2 2:2
chr10:-:19714521 0:1 0:0 0:0 0:0 2:2 0:2 0:2
chr10:-:19722422 1:1 0:0 0:0 0:0 2:2 2:2 0:0
chr10:-:19722588 0:1 0:0 0:0 0:0 0:2 0:2 0:0
```
.site deposits the phenotype of each splice site with both junction reads supporting the usage of the splice site and total junction reads in each sample.

After obtaining initial phenotype for each site (.site), we perform phenotype filtering 
3)
```
splice_read_count=*.site  (the initial phenotype file obtained from previous step)

out_file=*.filtered (deposit filtered phenotype file)
prop_file=*.prop
sd=0.01 (splice sites with variations less than the threshold will be deleted)
na_prop=0.4 (splice sites lacking phenotypes in larger than the threshold of the whole samples will be deleted)

$ISSAC pheno_output $splice_read_count $out_file $prop_file $sd $na_prop
```

The output of filter including two files:
(1).filtered (2).prop

## Sites located in single-intron clusters
For those sites without competing junction reads, we could quantify nonsplit reads crossing the splice sites to check if intron retention events exist around the splice site.  
ISSAC IR part provides the option for you to obtain nonsplit reads crossing the sites you are interested in.
```
bamfile=*.bam ###need original bam file with its corresponding index file
barcode=*.barcode  ###need the cell barcode of the cells you wanted to include
site=*.site  ###deposit the splice sites you are interested in
output_file=*.nonsplit  ###output CB-UMI based nonsplit reads crossing the splice site

$ISSAC IR extract -s FR -b ${barcode} -t ${site} ${bamfile} -o ${output_file}
```

## model construction
After obtaining phenotype(.filtered), GRM, genotype(.bcf &.csi),  .PC file, we could perform model construction for each splice site to estimate coefficients of fixed effect and random effect
GRM preparation

use plink to obtain .grm.bin, .grm.grm.N.bin file 
```
plink --bfile pruned_output --make-grm-bin --out grm_output
```

Here, we provided one approach to transform *.bin file to *.txt file for our downstream tasks.
```
library(plinkFile)

grm_bin<-"/home/e0950183/project/Imputation_April2024/filter/pruned/grm.grm.bin"
grm_N<-"/home/e0950183/project/Imputation_April2024/filter/pruned/grm.grm.N.bin"



prefix<-"/home/e0950183/project/Imputation_April2024/filter/pruned/grm"

dat<-readGRM(prefix)

newdat<-as.data.frame(dat)

write.table(newdat,"/home/e0950183/project/AIDA_phase2_cell/GRM/GRM.txt",sep=" ",row.names=TRUE,col.names=TRUE,quote=FALSE)
```
The final GRM file is shown as below:
```
JP_RIK_H001 JP_RIK_H002 JP_RIK_H003 JP_RIK_H004 JP_RIK_H005 JP_RIK_H006
JP_RIK_H001 0.989882528781891 0.0309822633862495 0.0303280465304852 0.0216292757540941 0.0148566737771034
JP_RIK_H002 0.0309822633862495 0.982672274112701 0.0118743935599923 0.0214019250124693 0.0196103043854237
JP_RIK_H003 0.0303280465304852 0.0118743935599923 0.986594617366791 0.022731801494956 0.00836068112403154
JP_RIK_H004 0.0216292757540941 0.0214019250124693 0.022731801494956 0.981386125087738 0.0213668830692768
JP_RIK_H005 0.0148566737771034 0.0196103043854237 0.00836068112403154 0.0213668830692768 0.982053458690643
JP_RIK_H006 0.0188088770955801 0.009139952249825 0.0207488182932138 0.0201078187674284 0.0271863844245672
```

Examples of PC file: Note: all the phenotypes and PCs’ names should be composed of ${sample}:${seq}; and ${sample} must exist in GRM file
```
JP_RIK_H002:1	JP_RIK_H002:3	JP_RIK_H003:1	JP_RIK_H003:10	JP_RIK_H003:14	JP_RIK_H003:3
-2.78345826143274	-2.54253101018537	-1.28028808817251	-1.65413653429238	-2.19730634684092	-2.03806144682479
1.36534270723467	0.625185320131184	-0.747500429552883	0.229720830781539	-0.0968962693380496	1.94665666355885
-1.64820219860581	-0.49179982238054	-0.382864839247662	-0.726471214189718	-0.823598633906727	-1.37327215340874
0.190717408902888	0.696397446346325	-0.787558003432512	-0.450713928140157	-1.06450187499565	-0.359833951076078
-0.246216596559398	-0.513259551637828	0.206514590906757	-0.0313668263851331	0.482904965754699	0.0864592524330276
-0.583672003889058	0.315994111242256	-0.3882240320452	0.190005776987925	0.429551644264448	0.636717879434517
```

ISSAC model part estimated parameters for both fixed effect and random effect through constructing null binomial mixed model.
```
$ISSAC model –s *.filtered   #phenotype    
       -p *.PC   # PC file  
       -g GRM.txt   # GRM file  
       -n 617 #  sample numbers of GRM file  
       -u ${model_pos} # the position of all the model files   
       -o  middle # the output name  
```

Note: ISSAC model will output all the splice site’s model containing in the phenotype file

Example of model files output by ISSAC model:

The first line includes: 1) splice site 2) dispersion parameter 3) variance components of random effects
another four lines are 1) residuals-null 2) π-null 3) total CB-UMI counts for the site in each sample 4)  CB-UMI counts supporting the usage of the site in each sample

## score tests for sQTL mapping

After obtaining model files, then we could utilize genotype file (*.bcf file with *.csi index file) to perform cis-sQTL mapping

ISSAC QTL mapping:
```
$ISSAC QTL   
       -p $model_file_pos  ## the model files position  
       -m  *.common ## the sample names of phenotype and PC  file  
       -c chr${i} ## the chromosome of splice sites  
       -x *.PC   ##PC file  
       -v chr${i}.recode.bcf ## genotype file  
       -w 1000000  ##SNPs located within the window will be used to perform cis-sQTL mapping  
       -o $output_pos   ##the output file’s position  
       -t 1  ##the sQTLs with pvalue less than this threshold will be output to result files  
```

Output of ISSAC QTL includes .result file and it will output results of all the splice sites in your desired chromosome of your .site file

Example of result file: 1) splice site 2) SNP 3) pvalue 4) effect size 5) standard error of effect size

```
```










