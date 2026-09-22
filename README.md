# ISSAC — Integrative single-cell splicing analysis and QTL caller

ISSAC is an integrated toolkit for single-cell splicing quantitative trait loci (sQTL) mapping and differential splicing analysis.

---

## Overview
ISSAC provides a complete pipeline for single-cell splicing analysis, including:
- Single-cell level junction extraction
- Site-based splice event quantification
- Null binomial model construction
- Score tests for cis-sQTL mapping

---

## Installation & Test Data

ISSAC can be run in three ways: natively on Linux (after building from source), via Docker (for Mac, Windows, and Linux), or via Apptainer/Singularity (for HPC clusters). A small example dataset is provided under `test_data/` to verify your installation and explore ISSAC's functionality.

### Clone the repository

```bash
git clone -b master https://github.com/boxiangliulab/ISSAC.git
cd ISSAC/test_data/
```

### Option 1: Native Linux build

#### Dependencies

The following C++ libraries must be available in your environment before building ISSAC (python version:3.8) if you want to compile yourself:

1. htslib 1.3
2. gsl
3. eigen3
4. nlopt
5. crypto (openssl)

```bash
conda install -c conda-forge gsl eigen nlopt openssl
conda install -c bioconda htslib=1.3
```
If you have built ISSAC from source on Linux (see the build instructions above), run the test pipeline directly:

```bash
bash linux_test.sh
```

### Option 2: Docker (Mac, Windows, Linux)

No local build required — this pulls a prebuilt container supporting both `amd64` and `arm64` (Apple Silicon) architectures.

```bash
bash docker_test.sh
```

> **Windows users:** run this script inside Git Bash or WSL, not PowerShell/CMD.

A successful run will show intermediate output such as `Converged`, `Genotype read in successfully`, and the pipeline completing through `Start pvalue computing!` without errors.

If you prefer to run ISSAC manually rather than via the wrapper script:

```bash
docker pull yuntian1999/issac:v1.3
docker run --rm yuntian1999/issac:v1.3 -h
```

To use your own data, mount your working directory and set it as the container's working directory:

```bash
docker run --rm -v "$PWD":"$PWD" -w "$PWD" yuntian1999/issac:v1.3 QTL [options]
```

### Option 3: Apptainer/Singularity (HPC users)

```bash
module load apptainer
bash apptainer_test.sh
```

If you prefer to run ISSAC manually:

```bash
apptainer pull issac.sif docker://yuntian1999/issac:v1.3
apptainer exec issac.sif ISSAC -h
```

To use your own data, bind mount your working directory:

```bash
apptainer exec --bind "$PWD":"$PWD" --pwd "$PWD" issac.sif ISSAC QTL [options]
```

---

## Notes

- Docker images are published at [`yuntian1999/issac`](https://hub.docker.com/r/yuntian1999/issac) and support both `linux/amd64` and `linux/arm64` platforms.
- The Apptainer `.sif` container is built directly from the same Docker image, so results are identical across all three installation methods.


### What the test pipeline covers
- Single-cell junction extraction
- Site-based splice event quantification
- Null binomial model construction
- cis-sQTL mapping

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
| `-o` | Output junction files |

### 1b. Per-Metacell Junction Statistics

Aggregates UMI-collapsed junction counts for each metacell using its associated barcode list.

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
| `-o` | Output statistics file (Barcode-UMI-unique read counts per junction per metacell) |

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
| `-t` | List of splice sites to quantify non-split read coverage |
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

Combines per-metacell non-split read files with junction read counts to generate site-level phenotypes for intron-retention-related splicing events. The output records splice-supporting reads and total reads (splice-supporting + non-split reads) at each site.

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
| `-l` | Splice-site list file |
| `-i` | File for total intron read counts across all metacells |
| `-o` | Output prefix for IR phenotype files |

### 2c. Phenotype Filtering

Filters out splice sites with low variability or high sparsity, retaining only informative sites for QTL mapping. Apply to both competitive intron and IR phenotype files. Splice-site usage variability is calculated using non-missing observations. For retained sites, metacells with zero total read count are imputed using the median numerator and denominator counts across non-missing metacells.

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
| `-s` | Minimum standard deviation threshold for splice-site usage |
| `-n` | Maximum missing-data fraction (sites with more zero total read count are treated as missing) |

---

## Step 3: Model Construction & QTL Mapping

Fits a binomial mixed model (GLMM) per splice site using a genetic relatedness matrix (GRM) to account for genetic relatedness and population structure, then performs cis-sQTL mapping within a defined window around each site.

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
  -t 10 \
  -v 0.05 \
  -i 30 \
  -l 0.001
```

| Flag | Description                                                                           |
| ---- | ------------------------------------------------------------------------------------- |
| `-s` | Filtered splicing phenotype file                                                      |
| `-p` | Principal components (PCs) file for covariate correction                              |
| `-n` | Number of individuals in the GRM file                                                 |
| `-g` | Genetic relatedness matrix (GRM) file for modeling genetic relatedness                |
| `-u` | Output directory/prefix for fitted null model files                                   |
| `-t` | Number of normalization parameter estimation iterations (×100)                        |
| `-v` | GRM sparsification threshold; relatedness values below this threshold are set to zero |
| `-i` | Maximum number of iterations for estimating fixed and random effects                  |
| `-l` | Convergence threshold for iterative parameter estimation (default: 0.001)             |


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
