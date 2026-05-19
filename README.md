# ISSAC — Integrative single-cell splicing analysis and QTL caller

ISSAC is an integrated toolkit for single-cell splicing quantitative trait loci (sQTL) mapping and differential splicing analysis.

---

## Overview
ISSAC provides a complete pipeline for single-cell splicing analysis, including:
- Single-cell level junction extraction
- Site-based splice event quantification
- Null binomial model construction
- Score tests for cis- and trans-sQTL mapping
- Differential splicing analysis

---

## Test Data

A small example dataset is provided to verify your installation and explore ISSAC's functionality.

### Download & Run
```bash
# Clone the repository
git clone -b master https://github.com/boxiangliulab/ISSAC.git

# Navigate to the test data directory
cd ISSAC/test_data/

# Run the test pipeline
bash test_detail.sh
```

### What the test pipeline covers
- Single-cell junction extraction
- Site-based splice event quantification
- Null binomial model construction
- cis- and trans-sQTL mapping
- Differential splicing analysis

### Expected output
```
test_data/
├── junctions_nonsplit_extract/          # Single-cell junction extraction
├── splice_phenotype_prepare/         # Site-based splice event quantification
├── model_construct_QTL_mapping/         # Null binomial model construction & downstream analysis
└── test_detail.sh  # Test script
```
---

## Table of Contents

- [Installation](#installation)
- [Pipeline Overview](#pipeline-overview)
- [Step 1: Junction & Non-split Read Extraction](#step-1-junction--non-split-read-extraction)
- [Step 2: Phenotype Preparation](#step-2-phenotype-preparation)
- [Step 3: Model Construction & QTL Mapping](#step-3-model-construction--qtl-mapping)

---

## Installation

```bash
git clone https://github.com/your-org/ISSAC.git
cd ISSAC/build
./ISSAC -h
```

Set the path to the compiled binary:

```bash
ISSAC=../build/ISSAC
```

---

## Pipeline Overview

```
BAM file
   │
   ▼
Step 1: Junction & non-split read extraction
   │
   ▼
Step 2: Phenotype preparation & filtering
   │
   ▼
Step 3: Null model fitting → cis-sQTL mapping
```

---

## Step 1: Junction & Non-split Read Extraction

Extract splicing junctions and non-split reads from BAM files. These provide raw per-barcode evidence for splicing events and intron retention.

### 1a. Junction Extraction

Extracts splice junctions from a 10X 3′ RNA-seq BAM file.

```bash
bamfile=junctions_nonsplit_extract/test.bam
output_junc=junctions_nonsplit_extract/test.junc

$ISSAC junctools extract \
  -a 8 \
  -m 50 \
  -M 500000 \
  -s RF \
  $bamfile \
  -o $output_junc
```

| Flag | Description |
|------|-------------|
| `-a` | Minimum anchor length (bp) |
| `-m` | Minimum intron length (bp) |
| `-M` | Maximum intron length (bp); junctions spanning longer distances are excluded |
| `-s` | Strand specificity (`RF` = reverse-forward, typical for 10X 3′ libraries) |

### 1b. Per-Metacell Junction Statistics

Aggregates junction read counts per barcode (metacell).

```bash
barcode=junctions_nonsplit_extract/metacell1.barcode
output_stat=junctions_nonsplit_extract/metacell1.stat

$ISSAC juncstat \
  -b $barcode \
  -j $output_junc \
  -o $output_stat
```

| Flag | Description |
|------|-------------|
| `-b` | Barcode list file defining which barcodes belong to this metacell |
| `-j` | Junction file produced in the previous step |
| `-o` | Output statistics file (read counts per junction per metacell) |

### 1c. Non-split Read Extraction

Extracts reads that do **not** span a splice junction, used for intron retention (IR) quantification.

```bash
site_list=junctions_nonsplit_extract/test.site
output_nonsplit=junctions_nonsplit_extract/metacell1.nonsplit

$ISSAC IR extract \
  -s RF \
  -b $barcode \
  -t $site_list \
  -a $bamfile \
  -o $output_nonsplit
```

| Flag | Description |
|------|-------------|
| `-s` | Strand specificity |
| `-b` | Barcode list for the metacell |
| `-t` | List of intronic sites to quantify non-split read coverage |
| `-a` | Input BAM file |
| `-o` | Output non-split read file |

---

## Step 2: Phenotype Preparation

Constructs per-metacell splicing phenotype matrices from junction and non-split read files. Two types of splicing events are handled:

- **(a) Competitive introns** — splice sites within intron clusters where multiple sites compete (classic sQTL signal)
- **(b) Single intron clusters** — intron retention events quantified from non-split reads

### 2a. Competitive Intron Phenotype Preparation

Groups splice sites across metacells and generates splice site usage ratios.

```bash
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
```

| Flag | Description |
|------|-------------|
| `-s` | File listing sample (metacell) names |
| `-j` | Directory containing per-metacell `.stat` files |
| `-o` | Output prefix for phenotype files |
| `-t` | Minimum total reads per junction across metacells |
| `-l` | Log file |
| `-n` | Minimum intron length (bp); should match junction extraction |
| `-x` | Maximum intron length (bp); should match junction extraction |

### 2b. Single Intron Cluster (Intron Retention) Phenotype Preparation

Combines per-metacell non-split read files into a site-level phenotype matrix, quantifying intron retention as the ratio of non-split to total reads at each site.

```bash
single_intron_site=splice_phenotype_prepare/test_single_intron_site
total_intron=splice_phenotype_prepare/phenotype_file/test.intron.out

$ISSAC IR_combine \
  -s $sample \
  -f splice_phenotype_prepare/nonsplit_file \
  -l $single_intron_site \
  -i $total_intron \
  -o splice_phenotype_prepare/phenotype_file/test
```

| Flag | Description |
|------|-------------|
| `-s` | Sample list file |
| `-f` | Directory containing per-metacell non-split read files |
| `-l` | List of intronic sites to include |
| `-i` | File for total intron read counts across all metacells |
| `-o` | Output prefix for IR phenotype files |

### 2c. Phenotype Filtering

Filters out splice sites with low variance or high sparsity, retaining only informative sites for QTL mapping. Apply to both competitive intron and IR phenotype files.

**Competitive intron sites:**

```bash
s=0.1
n=0.5

$ISSAC pheno_output \
  -r splice_phenotype_prepare/phenotype_file/test.site \
  -o splice_phenotype_prepare/phenotype_file/test.filtered \
  -p splice_phenotype_prepare/phenotype_file/test.prop \
  -s $s \
  -n $n
```

**Intron retention sites:**

```bash
$ISSAC pheno_output \
  -r splice_phenotype_prepare/phenotype_file/test_IR.site \
  -o splice_phenotype_prepare/phenotype_file/test_IR.filtered \
  -p splice_phenotype_prepare/phenotype_file/test_IR.prop \
  -s $s \
  -n $n
```

| Flag | Description |
|------|-------------|
| `-r` | Input site phenotype file |
| `-o` | Filtered output phenotype file |
| `-p` | Output file for per-site usage proportions |
| `-s` | Minimum variance threshold (sites below this are excluded) |
| `-n` | Maximum sparsity threshold (sites with more zero/missing values than this fraction are excluded) |

---

## Step 3: Model Construction & QTL Mapping

Fits a binomial mixed model (GLMM) per splice site using a genetic relatedness matrix (GRM) to control for population structure, then performs cis-sQTL mapping within a defined window around each site.

### 3a. Null Model Construction

Pre-fits null GLMMs (without genotype) for each splice site to avoid redundant computation during QTL mapping.

```bash
site_pheno=model_construct_QTL_mapping/gdT_GZMBhi_meta5_test.filtered
PC_file=model_construct_QTL_mapping/gdT_GZMBhi_meta5.PC
common_name=model_construct_QTL_mapping/gdT_GZMBhi_meta5.common

$ISSAC model \
  -s $site_pheno \
  -p $PC_file \
  -n 617 \
  -g model_construct_QTL_mapping/GRM.txt \
  -u model_construct_QTL_mapping/model \
  -t 10
```

| Flag | Description |
|------|-------------|
| `-s` | Filtered splicing phenotype file |
| `-p` | Genotype principal components (PCs) file for covariate correction |
| `-n` | Number of individuals in the GRM file |
| `-g` | Genetic relatedness matrix (GRM) file for controlling population stratification |
| `-u` | Output directory/prefix for fitted null model files |
| `-t` | Number of normalization parameter estimation iterations (×10) |

Collect sites for which null models were successfully built:

```bash
ls model_construct_QTL_mapping/model/* \
  | cut -d '/' -f 3 \
  | cut -d '.' -f 1 \
  > model_construct_QTL_mapping/test_site.list
```

### 3b. cis-sQTL Mapping

Tests association between each splice site and all SNPs within a cis window, using pre-fitted null models for efficient GLMM-based testing. Process one chromosome at a time to enable parallelisation.

```bash
genotype=model_construct_QTL_mapping/test.bcf
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
```

| Flag | Description |
|------|-------------|
| `-s` | List of splice sites to test |
| `-o` | Output prefix for QTL result files |
| `-c` | Chromosome to map |
| `-v` | Genotype file in BCF format (must be indexed) |
| `-x` | PC file (same covariates used during model construction) |
| `-p` | Directory containing pre-fitted null model files |
| `-w` | cis window size in bp (±500 kb around each splice site) |
| `-m` | File listing sample names in the same order as the null model files |
| `-t` | P-value threshold; associations with p > threshold are discarded |
> **Note:** Please ensure all dependencies are installed before running the test. Refer to the [Installation](#installation) section for details.

---

## License

---

## Contact
- **Lab**: [Boxiang Liu Lab](https://github.com/boxiangliulab)
- **Issues**: Please report bugs via the [GitHub Issues](https://github.com/boxiangliulab/ISSAC/issues) page
