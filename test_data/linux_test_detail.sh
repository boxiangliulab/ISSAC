###############################################################
# ISSAC Test Pipeline: Integrative Single-cell Splicing Analysis and QTL Caller
# This script demonstrates the ISSAC workflow using test data,
# covering junction extraction, phenotype preparation, null-model
# construction, and cis-sQTL mapping.
###############################################################

ISSAC=../build/ISSAC


###############################################################
# Step 1: Junction & Non-split Read Extraction
# Extract splice junctions and non-split reads from BAM files.
# Junction reads retain cell-barcode and UMI information for
# molecule-level deduplication.
###############################################################

# --- Junction Extraction (10X 3' RNA-seq BAM) ---
# Extract splice junctions from a BAM file.
# -a 8        : Minimum anchor length (bp)
# -m 50       : Minimum intron length (bp); shorter junctions are excluded
# -M 500000   : Maximum intron length (bp); longer junctions are excluded
# -s RF       : Strand specificity (RF = reverse-forward, typical for 10X 3' libraries)

bamfile=junctions_nonsplit_extract/test.bam
output_junc=junctions_nonsplit_extract/test.junc

$ISSAC junctools extract \
  -a 8 \
  -m 50 \
  -M 500000 \
  -s RF \
  $bamfile \
  -o $output_junc


# --- Per-Metacell Junction Statistics ---
# Aggregates UMI-collapsed junction counts for a metacell using
# its associated barcode list.
# -b : Barcode list file defining which barcodes belong to the metacell
# -j : Junction file produced in the previous step
# -o : Output statistics file containing UMI-collapsed junction counts

barcode=junctions_nonsplit_extract/metacell1.barcode
output_stat=junctions_nonsplit_extract/metacell1.stat

$ISSAC juncstat \
  -b $barcode \
  -j $output_junc \
  -o $output_stat


# --- Non-split Read Extraction ---
# Extract reads crossing specified splice sites without splicing.
# These reads are used for intron-retention-related quantification.
# -s : Strand specificity
# -b : Barcode list for the metacell
# -t : List of splice sites to quantify non-split read coverage
# -a : Input BAM file
# -o : Output non-split read file

site_list=junctions_nonsplit_extract/test.site
output_nonsplit=junctions_nonsplit_extract/metacell1.nonsplit

$ISSAC IR extract \
  -s RF \
  -b $barcode \
  -t $site_list \
  -a $bamfile \
  -o $output_nonsplit


###############################################################
# Step 2: Phenotype Preparation
# Construct per-metacell splicing phenotype matrices from junction
# and non-split read files. Two types of splicing events are handled:
#   (a) Competitive intron clusters
#   (b) Intron-retention-related events incorporating non-split reads
###############################################################

# --- (a) Competitive Intron Phenotype Preparation ---
# Group junctions into intron clusters and quantify splice-site usage.
# -s : File listing sample (metacell) names
# -j : Directory containing per-metacell .stat files
# -o : Output prefix for phenotype files
# -t : Minimum total junction read count across all metacells
# -l : Log file
# -n : Minimum intron length (bp); consistent with junction extraction
# -x : Maximum intron length (bp); consistent with junction extraction

ls splice_phenotype_prepare/stat_file/*.stat \
  | cut -d '/' -f 3 \
  | cut -d '.' -f 1 \
  > splice_phenotype_prepare/sample_file

sample=splice_phenotype_prepare/sample_file

$ISSAC pheno_group \
  -s $sample \
  -j splice_phenotype_prepare/stat_file \
  -o splice_phenotype_prepare/phenotype_file/test \
  -t 50 \
  -l log.out \
  -n 50 \
  -x 500000


# --- (b) Intron-Retention-Related Phenotype Preparation ---
# Combine junction counts with per-metacell non-split read counts
# to generate site-level phenotypes. The output records splice-supporting
# reads as the numerator and total reads (splice-supporting + non-split)
# as the denominator.
# -s : Sample list file
# -f : Directory containing per-metacell non-split read files
# -l : Splice-site list file
# -i : Intron read-count file generated above
# -o : Output prefix for IR phenotype files

single_intron_site=splice_phenotype_prepare/test_single_intron_site
total_intron=splice_phenotype_prepare/phenotype_file/test.intron.out

$ISSAC IR_combine \
  -s $sample \
  -f splice_phenotype_prepare/nonsplit_file \
  -l $single_intron_site \
  -i $total_intron \
  -o splice_phenotype_prepare/phenotype_file/test


# --- Phenotype Filtering: Competitive Intron Sites ---
# Filter splice sites with low usage variability or high missingness.
# Usage standard deviation is calculated using non-missing observations.
# For retained sites, metacells with zero total read count are imputed
# using the median numerator and denominator counts across non-missing
# metacells.
#
# -r     : Input site phenotype file
# -o     : Filtered output phenotype file
# -p     : Output file for per-site usage proportions
# -s 0.1 : Minimum standard deviation threshold for splice-site usage
# -n 0.5 : Maximum missing-data fraction; metacells with zero total
#          read count are treated as missing

site_file=splice_phenotype_prepare/phenotype_file/test.site
filter_output=splice_phenotype_prepare/phenotype_file/test.filtered
filter_prop=splice_phenotype_prepare/phenotype_file/test.prop

s=0.1
n=0.5

$ISSAC pheno_output \
  -r $site_file \
  -o $filter_output \
  -p $filter_prop \
  -s $s \
  -n $n


# --- Phenotype Filtering: Intron-Retention-Related Sites ---
# Apply the same filtering and imputation procedure.

site_file=splice_phenotype_prepare/phenotype_file/test_IR.site
filter_output=splice_phenotype_prepare/phenotype_file/test_IR.filtered
filter_prop=splice_phenotype_prepare/phenotype_file/test_IR.prop

$ISSAC pheno_output \
  -r $site_file \
  -o $filter_output \
  -p $filter_prop \
  -s $s \
  -n $n


###############################################################
# Step 3: Model Construction & QTL Mapping
# Fit a binomial mixed model (GLMM) for each splice site using a
# genetic relatedness matrix (GRM) to account for genetic relatedness
# and population structure, followed by cis-sQTL mapping.
###############################################################

site_pheno=model_construct_QTL_mapping/gdT_GZMBhi_meta5_test.filtered
PC_file=model_construct_QTL_mapping/gdT_GZMBhi_meta5.PC
genotype=model_construct_QTL_mapping/test.bcf
common_name=model_construct_QTL_mapping/gdT_GZMBhi_meta5.common


# --- Model Construction ---
# Fit null GLMMs (without genotype) for each splice site.
# Pre-fitting null models avoids redundant computation during QTL mapping.
#
# -s       : Filtered splicing phenotype file
# -p       : Principal components (PCs) file for covariate correction
# -n 617   : Number of individuals in the GRM file
# -g       : Genetic relatedness matrix (GRM) file
# -u       : Output directory/prefix for fitted null model files
# -t 10    : Number of normalization parameter estimation iterations (x100)
# -v 0.05  : GRM sparsification threshold; relatedness values below
#            this threshold are set to zero
# -i 30    : Maximum number of iterations for estimating fixed and
#            random effects
# -l 0.001 : Convergence threshold for iterative parameter estimation

$ISSAC model \
  -s $site_pheno \
  -p $PC_file \
  -n 617 \
  -g model_construct_QTL_mapping/GRM.txt \
  -u model_construct_QTL_mapping/model \
  -t 10 \
  -v 0.05 \
  -i 30 \
  -l 0.001


# Collect sites for which null models were successfully built

ls model_construct_QTL_mapping/model/* \
  | cut -d '/' -f 3 \
  | cut -d '.' -f 1 \
  > model_construct_QTL_mapping/test_site.list


# --- cis-sQTL Mapping ---
# Test association between each splice site and SNPs within a cis window
# using the pre-fitted null models.
#
# -s        : List of splice sites to test
# -o        : Output directory/prefix for QTL result files
# -c        : Chromosome to map
# -v        : Genotype file in BCF format (must be indexed)
# -x        : PC file used during model construction
# -p        : Directory containing pre-fitted null model files
# -w 500000 : cis window size (±500 kb around each splice site)
# -m        : File listing sample names in the same order as the null model
# -t        : P-value threshold; associations with P > threshold are discarded

site_list=model_construct_QTL_mapping/test_site.list
chr=chr10

$ISSAC QTL \
  -s $site_list \
  -o model_construct_QTL_mapping/qtl \
  -c $chr \
  -v $genotype \
  -x $PC_file \
  -p model_construct_QTL_mapping/model \
  -w 500000 \
  -m $common_name \
  -t 1

