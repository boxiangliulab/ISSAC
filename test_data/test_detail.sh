###############################################################
# ISSAC Test Pipeline: Integrative Single-cell Splicing Analysis and qtl Caller
# This script demonstrates the full ISSAC workflow using test data,
# covering junction extraction, phenotype preparation, and QTL mapping.
###############################################################

ISSAC=../build/ISSAC

###############################################################
# Step 1: Junction & Non-split Read Extraction
# Extract splicing junctions and non-split reads from BAM files.
# These form the raw evidence for splicing events per cell barcode.
###############################################################

# --- Junction Extraction (10X 3' RNA-seq BAM) ---
# Extracts splice junctions from a BAM file.
# -a 8        : Minimum anchor length (bp)
# -m 50       : Maximum intron length (bp); junctions spanning shorter distances are excluded
# -M 500000   : Maximum intron length (bp); junctions spanning longer distances are excluded
# -s RF       : Strand specificity (RF = reverse-forward, typical for 10X 3' libraries)
bamfile=junctions_nonsplit_extract/test.bam
output_junc=junctions_nonsplit_extract/test.junc

$ISSAC junctools extract -a 8 -m 50 -M 500000 -s RF $bamfile -o $output_junc

# --- Per-Metacell Junction Statistics ---
# Aggregates junction read counts per barcode (metacell).
# -b : Barcode list file defining which barcodes belong to this metacell
# -j : Junction file produced in the previous step
# -o : Output statistics file (read counts per junction per metacell)
barcode=junctions_nonsplit_extract/metacell1.barcode
output_stat=junctions_nonsplit_extract/metacell1.stat

$ISSAC juncstat -b $barcode -j $output_junc -o $output_stat

# --- Non-split Read Extraction ---
# Extracts reads that do NOT span a splice junction (used for intron retention quantification).
# -s     : Strand specificity
# -b     : Barcode list for the metacell
# -t     : List of intronic sites to quantify non-split read coverage
# -a     : Input BAM file
# -o     : Output non-split read file
site_list=junctions_nonsplit_extract/test.site
output_nonsplit=junctions_nonsplit_extract/metacell1.nonsplit

$ISSAC IR extract -s RF -b $barcode -t $site_list -a $bamfile -o $output_nonsplit


###############################################################
# Step 2: Phenotype Preparation
# Constructs per-metacell splicing phenotype matrices from junction
# and non-split read files. Two types of splicing events are handled:
#   (a) Competitive introns: splice sites within intron clusters
#       where multiple splice sites compete (classic sQTL signal)
#   (b) Single intron clusters: intron retention events quantified
#       from non-split reads
###############################################################

# --- (a) Competitive Intron Phenotype Preparation ---
# Groups splice sites into intron clusters across metacells.
# Generates splice site usage ratios for competitive splice sites.
# -s : File listing sample (metacell) names
# -j : Directory containing per-metacell .stat files
# -o : Output prefix for phenotype files
# -t : Minimum total reads per junction across all metacells
# -l : Log file
# -n : Minimum intron length (bp); consistent with junction extraction
# -x : Maximum intron length (bp); consistent with junction extraction
ls splice_phenotype_prepare/stat_file/*.stat | cut -d '/' -f 3 | cut -d '.' -f 1 \
  > splice_phenotype_prepare/sample_file
sample=splice_phenotype_prepare/sample_file

$ISSAC pheno_group -s $sample \
  -j splice_phenotype_prepare/stat_file \
  -o splice_phenotype_prepare/phenotype_file/test \
  -t 50 \
  -l log.out \
  -n 50 \
  -x 500000

# --- (b) Single Intron Cluster (Intron Retention) Phenotype Preparation ---
# Combines per-metacell non-split read files into a site-level phenotype matrix.
# Quantifies intron retention as the ratio of split reads to total reads (including split and non-split reads) at each site.
# -s : Sample list file
# -f : Directory containing per-metacell non-split read files
# -l : List of intronic sites to include
# -i : File for total intron read counts across all metacells
# -o : Output prefix for IR phenotype files
single_intron_site=splice_phenotype_prepare/test_single_intron_site
total_intron=splice_phenotype_prepare/phenotype_file/test.intron.out

$ISSAC IR_combine -s $sample \
  -f splice_phenotype_prepare/nonsplit_file \
  -l $single_intron_site \
  -i $total_intron \
  -o splice_phenotype_prepare/phenotype_file/test

# --- Phenotype Filtering: Competitive Intron Sites ---
# Filters out splice sites with low variance or high sparsity across metacells,
# retaining only informative sites for QTL mapping.
# -r : Input site phenotype file
# -o : Filtered output phenotype file
# -p : Output file for per-site usage proportions
# -s 0.1 : Minimum variance threshold (sites with variance < 0.1 are excluded)
# -n 0.5 : Maximum sparsity threshold (sites with >50% zero/missing values are excluded)
site_file=splice_phenotype_prepare/phenotype_file/test.site
filter_output=splice_phenotype_prepare/phenotype_file/test.filtered
filter_prop=splice_phenotype_prepare/phenotype_file/test.prop
s=0.1
n=0.5

$ISSAC pheno_output -r $site_file \
  -o $filter_output \
  -p $filter_prop \
  -s $s \
  -n $n

# --- Phenotype Filtering: Single Intron Cluster (IR) Sites ---
# Same filtering logic applied to the intron retention phenotype matrix.
site_file=splice_phenotype_prepare/phenotype_file/test_IR.site
filter_output=splice_phenotype_prepare/phenotype_file/test_IR.filtered
filter_prop=splice_phenotype_prepare/phenotype_file/test_IR.prop

$ISSAC pheno_output -r $site_file \
  -o $filter_output \
  -p $filter_prop \
  -s $s \
  -n $n


###############################################################
# Step 3: Model Construction & QTL Mapping
# Fits a binomial mixed model (GLMM) for each splice site using a
# genetic relatedness matrix (GRM) to control for population structure,
# then performs cis-sQTL mapping within a defined window around each site.
###############################################################

site_pheno=model_construct_QTL_mapping/gdT_GZMBhi_meta5_test.filtered
PC_file=model_construct_QTL_mapping/gdT_GZMBhi_meta5.PC
genotype=model_construct_QTL_mapping/test.bcf
common_name=model_construct_QTL_mapping/gdT_GZMBhi_meta5.common

# --- Model Construction ---
# Fits null GLMMs (without genotype) for each splice site.
# Pre-fitting null models avoids redundant computation during QTL mapping.
# -s : Filtered splicing phenotype file
# -p : Genotype principal components (PCs) file for covariate correction
# -n : Number of individuals in GRM file
# -g : Genetic relatedness matrix (GRM) file for controlling population stratification
# -u : Output directory/prefix for fitted null model files
# -t : Times of normalization parameters estimation (×10)
$ISSAC model -s $site_pheno \
  -p $PC_file \
  -n 617 \
  -g model_construct_QTL_mapping/GRM.txt \
  -u model_construct_QTL_mapping/model \
  -t 10

# Collect the list of sites for which null models were successfully built
ls model_construct_QTL_mapping/model/* | cut -d '/' -f 3 | cut -d '.' -f 1 \
  > model_construct_QTL_mapping/test_site.list

# --- cis-sQTL Mapping ---
# Tests association between each splice site and all SNPs within a cis window.
# Uses the pre-fitted null models for efficient GLMM-based association testing.
# -s : List of splice sites to test
# -o : Output directory for QTL result files
# -c : Chromosome to map (process one chromosome at a time for parallelisation)
# -v : Genotype file in BCF format (must be indexed)
# -x : PC file (same covariates used during model construction)
# -p : Directory containing pre-fitted null model files
# -w 500000 : cis window size (±500 kb around each splice site)
# -m : File listing sample names (the same sequence as in null model files)
# -t : Threshold for significance (Associations with Pvalue > t will be discarded)
site_list=model_construct_QTL_mapping/test_site.list
chr=chr10

$ISSAC QTL -s $site_list \
  -o model_construct_QTL_mapping/qtl \
  -c $chr \
  -v $genotype \
  -x $PC_file \
  -p model_construct_QTL_mapping/model \
  -w 500000 \
  -m $common_name \
  -t 1
