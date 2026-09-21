### ISSAC test data (Apptainer version — for HPC users)

# Load Apptainer module (adjust to your HPC's module name if different)
module load apptainer

# Path to the ISSAC .sif container (pull once if not already present)
ISSAC_SIF=issac.sif
if [ ! -f "$ISSAC_SIF" ]; then
    apptainer pull "$ISSAC_SIF" docker://yuntian1999/issac:v1.3
fi

# Wrap `apptainer exec` as a function so the rest of the script can call
# ISSAC exactly like a local binary. We bind the current working
# directory into the container at the same path, so all relative paths
# below (junctions_nonsplit_extract/..., splice_phenotype_prepare/...,
# model_construct_QTL_mapping/...) work unchanged.
ISSAC() {
  apptainer exec --bind "$PWD":"$PWD" --pwd "$PWD" "$ISSAC_SIF" ISSAC "$@"
}


### Step 1: Junction & non-split read extraction

# Junction extraction from 10X 3' RNA-seq BAM
bamfile=junctions_nonsplit_extract/test.bam
output_junc=junctions_nonsplit_extract/test.junc

ISSAC junctools extract -a 8 -m 50 -M 500000 -s RF $bamfile -o $output_junc

# Junction statistics
barcode=junctions_nonsplit_extract/metacell1.barcode
output_stat=junctions_nonsplit_extract/metacell1.stat

ISSAC juncstat -b $barcode -j $output_junc -o $output_stat

# Non-split read extraction
site_list=junctions_nonsplit_extract/test.site
output_nonsplit=junctions_nonsplit_extract/metacell1.nonsplit

ISSAC IR extract -s RF -b $barcode -t $site_list -a $bamfile -o $output_nonsplit


### Step 2: Phenotype preparation & output
### Both competitive intron and single intron cluster types

## Prepare splice sites located in competitive introns
ls splice_phenotype_prepare/stat_file/*.stat \
  | cut -d '/' -f 3 \
  | cut -d '.' -f 1 \
  > splice_phenotype_prepare/sample_file

sample=splice_phenotype_prepare/sample_file

ISSAC pheno_group -s $sample \
  -j splice_phenotype_prepare/stat_file \
  -o splice_phenotype_prepare/phenotype_file/test \
  -t 50 \
  -l log.out \
  -n 50 \
  -x 500000

## Prepare intron-retention-related site phenotypes by combining
## junction and non-split read counts
single_intron_site=splice_phenotype_prepare/test_single_intron_site
total_intron_output=splice_phenotype_prepare/phenotype_file/test.intron.out

ISSAC IR_combine -s $sample \
  -f splice_phenotype_prepare/nonsplit_file \
  -l $single_intron_site \
  -i $total_intron_output \
  -o splice_phenotype_prepare/phenotype_file/test

## Filter sites in competitive introns with low variability and high sparsity
site_file=splice_phenotype_prepare/phenotype_file/test.site
filter_output=splice_phenotype_prepare/phenotype_file/test.filtered
filter_prop=splice_phenotype_prepare/phenotype_file/test.prop

s=0.1
n=0.5

ISSAC pheno_output -r $site_file \
  -o $filter_output \
  -p $filter_prop \
  -s $s \
  -n $n

## Filter sites in single intron clusters with low variability and high sparsity
site_file=splice_phenotype_prepare/phenotype_file/test_IR.site
filter_output=splice_phenotype_prepare/phenotype_file/test_IR.filtered
filter_prop=splice_phenotype_prepare/phenotype_file/test_IR.prop

ISSAC pheno_output -r $site_file \
  -o $filter_output \
  -p $filter_prop \
  -s $s \
  -n $n


### Step 3: Model construction & QTL mapping

site_pheno=model_construct_QTL_mapping/gdT_GZMBhi_meta5_test.filtered
PC_file=model_construct_QTL_mapping/gdT_GZMBhi_meta5.PC
genotype=model_construct_QTL_mapping/test.bcf
common_name=model_construct_QTL_mapping/gdT_GZMBhi_meta5.common

ISSAC model \
  -s $site_pheno \
  -p $PC_file \
  -n 617 \
  -g model_construct_QTL_mapping/GRM.txt \
  -v 0.05 \
  -i 30 \
  -l 0.001 \
  -u model_construct_QTL_mapping/model \
  -t 10

ls model_construct_QTL_mapping/model/* \
  | cut -d '/' -f 3 \
  | cut -d '.' -f 1 \
  > model_construct_QTL_mapping/test_site.list

site_list=model_construct_QTL_mapping/test_site.list
chr=chr10

ISSAC QTL -s $site_list \
  -o model_construct_QTL_mapping/qtl \
  -c $chr \
  -v $genotype \
  -x $PC_file \
  -p model_construct_QTL_mapping/model \
  -w 500000 \
  -m $common_name \
  -t 1
