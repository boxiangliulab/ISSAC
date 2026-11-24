# ISSAC — Integrative single-cell splicing analysis and QTL caller

ISSAC is a toolkit for mapping splicing quantitative trait loci (sQTLs) in single-cell and single-nucleus RNA-seq data. It quantifies splice-site usage per metacell and tests genetic associations with a binomial generalized linear mixed model.

---

## Overview

ISSAC provides a complete pipeline for single-cell splicing analysis:

- Single-cell junction and non-split read extraction
- Site-based splice event quantification
- Null binomial mixed model construction
- Score tests for cis-sQTL mapping

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

## Table of Contents

- [Installation & Test Data](#installation--test-data)
- [Input conventions](#input-conventions)
- [Step 1: Junction & Non-split Read Extraction](#step-1-junction--non-split-read-extraction)
- [Step 2: Phenotype Preparation](#step-2-phenotype-preparation)
- [Step 3: Model Construction & QTL Mapping](#step-3-model-construction--qtl-mapping)
- [Contact](#contact)

---

## Installation & Test Data

ISSAC can be run in three ways: natively on Linux (after building from source), via Docker (Mac, Windows and Linux), or via Apptainer/Singularity (HPC clusters). A small example dataset is provided under `test_data/` to verify the installation.

### Option 1: Native Linux build

#### Requirements

- A C++ compiler (e.g. `g++`) and CMake (>= 3.10)
- htslib, gsl, eigen3, nlopt and openssl, installed via conda as shown below

#### Build and test

```bash
git clone -b master https://github.com/boxiangliulab/ISSAC.git
cd ISSAC

# Create and activate an environment with the dependencies
conda create -n issac -c conda-forge -c bioconda gsl eigen nlopt openssl htslib
conda activate issac

# Build (activate the environment first so that CMake can find the libraries)
rm -rf build && mkdir build && cd build
cmake .. && make

# Run the test pipeline
cd ../test_data && bash linux_test.sh
```

After building, set the path to the binary used in the examples below:

```bash
ISSAC=../build/ISSAC
$ISSAC -h
```

### Option 2: Docker (Mac, Windows, Linux)

No local build is required. The prebuilt image supports both `amd64` and `arm64` (Apple Silicon).

```bash
git clone -b master https://github.com/boxiangliulab/ISSAC.git
cd ISSAC/test_data
bash docker_test.sh
```

> **Windows users:** run this script in Git Bash or WSL, not PowerShell/CMD.

To run ISSAC manually:

```bash
docker pull yuntian1999/issac:v1.3
docker run --rm yuntian1999/issac:v1.3 -h
```

To use your own data, mount your working directory and set it as the container's working directory:

```bash
docker run --rm -v "$PWD":"$PWD" -w "$PWD" yuntian1999/issac:v1.3 QTL [options]
```

### Option 3: Apptainer/Singularity (HPC)

```bash
git clone -b master https://github.com/boxiangliulab/ISSAC.git
cd ISSAC/test_data
module load apptainer
bash apptainer_test.sh
```

To run ISSAC manually:

```bash
apptainer pull issac.sif docker://yuntian1999/issac:v1.3
apptainer exec issac.sif ISSAC -h
```

To use your own data, bind your working directory:

```bash
apptainer exec --bind "$PWD":"$PWD" --pwd "$PWD" issac.sif ISSAC QTL [options]
```

> On some HPC systems, user namespaces are disabled on login nodes; run Apptainer on a compute node.

### Notes

- Docker images are published at [`yuntian1999/issac`](https://hub.docker.com/r/yuntian1999/issac) for `linux/amd64` and `linux/arm64`.
- The Apptainer `.sif` container is built from the same Docker image, so results are identical across the three installation methods.

### What the test pipeline covers

A successful run prints intermediate messages such as `Converged`, `Genotype read in successfully` and `Start pvalue computing!` without errors, and produces:

```
test_data/
├── junctions_nonsplit_extract/      # Single-cell junction extraction
├── splice_phenotype_prepare/        # Site-based splice event quantification
├── model_construct_QTL_mapping/     # Null binomial model construction & cis-sQTL mapping
├── linux_test.sh                    # Test script for the native Linux build
├── docker_test.sh                   # Test script using Docker
└── apptainer_test.sh                # Test script using Apptainer/Singularity
```

---

## Input conventions

These conventions link metacells, donors, covariates and genotypes. Inputs that do not follow them may not raise an error but will give incorrect results.

**Metacell names.** Each metacell is named `<donor ID>:<index>`, e.g. `D001:0`, `D001:1`. The donor ID (the part before the first `:`) must match the donor IDs in the GRM file and the sample IDs in the genotype (VCF/BCF) file. Metacells from the same donor share the donor's genotype and random effect.

**GRM file (`-g`).** Space-delimited. The first line lists the donor IDs; each following line contains one donor ID followed by its relatedness values with all donors, in the same order as the first line:

```
D001 D002 D003
D001 1.002 0.003 -0.010
D002 0.003 0.998 0.001
D003 -0.010 0.001 1.005
```

The number of donors must equal the value passed to `-n`. The GRM may include donors that are absent from the phenotype data.

**Covariate (PC) file (`-p` in `model`, `-x` in `QTL`).** Tab-delimited. The first line lists metacell names; each following line contains one covariate (e.g. a principal component, sex or batch indicator) with one value per metacell and no row label. An intercept is added automatically. The same file must be used for model construction and QTL mapping.

```
D001:0	D001:1	D002:0	...
0.153	0.149	-0.021	...
1	1	0	...
```

Only metacells present in both the phenotype file and the covariate file are analyzed.

---

## Step 1: Junction & Non-split Read Extraction

Extract splicing junctions and non-split reads from BAM files. These provide per-barcode evidence for splicing events and intron retention.

### 1a. Junction extraction

Extracts splice junctions from a 10x 3′ RNA-seq BAM file.

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
| `-s` | Strand specificity (`RF` = reverse-forward, typical for 10x 3′ libraries) |
| `-o` | Output junction file |

### 1b. Per-metacell junction statistics

Aggregates UMI-collapsed junction counts for each metacell using its barcode list.

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
| `-b` | Barcode list defining the barcodes that belong to this metacell |
| `-j` | Junction file produced in step 1a |
| `-o` | Output statistics file (UMI-unique read counts per junction for this metacell) |

### 1c. Non-split read extraction

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
| `-t` | List of splice sites at which non-split read coverage is quantified |
| `-a` | Input BAM file |
| `-o` | Output non-split read file |

---

## Step 2: Phenotype Preparation

Constructs per-metacell splicing phenotype matrices from the junction and non-split read files. Two types of splicing events are handled:

- **(a) Competitive introns** — splice sites within intron clusters where multiple sites compete
- **(b) Single-intron clusters** — intron retention events quantified from non-split reads

### 2a. Competitive intron phenotype preparation

Groups splice sites across metacells and generates splice-site usage counts.

```bash
ls splice_phenotype_prepare/stat_file | sed -n 's/\.stat$//p' > splice_phenotype_prepare/sample_file
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
| `-s` | File listing metacell names, one per line |
| `-j` | Directory containing per-metacell `.stat` files |
| `-o` | Output prefix for phenotype files |
| `-t` | Minimum total reads per junction across metacells |
| `-l` | Log file |
| `-n` | Minimum intron length (bp); should match junction extraction |
| `-x` | Maximum intron length (bp); should match junction extraction |

### 2b. Single-intron cluster (intron retention) phenotype preparation

Combines per-metacell non-split read files with junction read counts to generate site-level phenotypes for intron-retention events. The output records splice-supporting reads and total reads (splice-supporting + non-split reads) at each site.

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
| `-s` | Metacell list file |
| `-f` | Directory containing per-metacell non-split read files |
| `-l` | Splice-site list file |
| `-i` | File of total intron read counts across all metacells |
| `-o` | Output prefix for IR phenotype files |

### 2c. Phenotype filtering

Removes splice sites with low variability or high sparsity, retaining informative sites for QTL mapping. Apply to both competitive intron and IR phenotype files. Splice-site usage variability is calculated from non-missing observations. For retained sites, metacells with zero total read count are imputed using the median numerator and denominator counts across non-missing metacells.

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
| `-p` | Output file of per-site usage proportions |
| `-s` | Minimum standard deviation of splice-site usage |
| `-n` | Maximum fraction of missing metacells (metacells with zero total read count are treated as missing) |

The filtered phenotype file is space-delimited: the first line lists metacell names, and each following line contains a splice site name followed by one `used:total` read count per metacell.

---

## Step 3: Model Construction & QTL Mapping

For each splice site, ISSAC models splice-site read counts with a binomial generalized linear mixed model (logit link). The model includes covariates, a donor-level random effect structured by the genetic relatedness matrix (GRM), which captures the correlation among metacells from the same donor as well as relatedness between donors, and an observation-level random effect for overdispersion. Null models are fitted once per site; cis-sQTL mapping then uses score tests against the pre-fitted null models.

### 3a. Null model construction

Fits the null GLMM (without genotype) for each splice site and estimates the site-specific variance correction factor *r*, which accounts for the discrepancy between the true variance of the score statistic and the variance under the naïve working-weight matrix. Two estimation modes are available:

| Mode | Option | Recommended for |
| ---- | ------ | --------------- |
| **Permutation** (default) | `-m permute` | Cohorts with limited between-donor relatedness (e.g., most population-based single-cell cohorts) |
| **Real-marker sampling** (SAIGE-style) | `-m sample -r <file>` | Cohorts with substantial between-donor relatedness |

#### Option 1: Permutation-based correction (default)

*r* is estimated by permuting a synthetic genotype vector (MAF = 0.5) across donors, with all metacells from the same donor sharing the same permuted genotype. This approach does not require stable estimation of the variance components and therefore applies to all splice sites, including those at which the mixed model fails to converge. It assumes that donors are approximately exchangeable under the null.

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

In this mode, splice sites are handled according to their model-fitting status:

| Model-fitting status | Model used for testing | Variance correction |
| -------------------- | ---------------------- | ------------------- |
| Mixed model converged | Binomial GLMM | Permutation-based *r* |
| Mixed model not converged | Fixed-effects binomial model (status `Fixed`) | Permutation-based *r*, which also absorbs the within-donor correlation not modelled by the fixed-effects model |
| Both models failed | — | Site excluded (no output) |

#### Option 2: SAIGE-style correction using real markers

*r* is estimated as the mean ratio G<sup>T</sup>PG / G<sup>T</sup>WG across user-provided real, unpermuted markers, which preserves the relatedness structure among donors.

```bash
$ISSAC model \
  -s $site_pheno \
  -p $PC_file \
  -n 617 \
  -g model_construct_QTL_mapping/GRM.txt \
  -u model_construct_QTL_mapping/model \
  -v 0.05 \
  -i 30 \
  -l 0.001 \
  -m sample \
  -r model_construct_QTL_mapping/null_markers.txt
```

In this mode, splice sites are handled according to their model-fitting status:

| Model-fitting status | Variance correction |
| -------------------- | ------------------- |
| Converged, τ<sub>g</sub> > 10<sup>-6</sup> | SAIGE-style G<sup>T</sup>PG / G<sup>T</sup>WG ratio |
| Not converged, or τ<sub>g</sub> at the lower boundary | Site excluded (no output) |

**Marker genotype file format (`-r`).** Tab-delimited; the first column is the donor ID, and each remaining column is one marker with genotype dosages (0/1/2; missing values as `NA`). An optional header line is allowed.

```
donor_id   rs1001   rs1002   rs1003   ...
D001       0        1        2
D002       1        0        1
D003       2        1        0
```

> **Important:** rows are matched to donors **by order**, not by the ID in the first column. Donors must be listed in the same order as in the GRM file.

We recommend using at least 30 common (MAF ≥ 0.05), approximately independent markers located on chromosomes other than those of the tested splice sites, so that they have no *cis* effect on splicing.

#### Options

| Flag | Description |
| ---- | ----------- |
| `-s` | Filtered splicing phenotype file |
| `-p` | Covariate (PC) file |
| `-n` | Number of donors in the GRM file |
| `-g` | Genetic relatedness matrix (GRM) file |
| `-u` | Output directory for fitted null model files |
| `-t` | Number of permutations for estimating *r*, in units of 100 (default: 10, i.e., 1,000 permutations); used only in permutation mode |
| `-v` | GRM sparsification threshold; relatedness values at or below this threshold are set to zero (recommended: 0.05) |
| `-i` | Maximum number of iterations for estimating fixed and random effects (recommended: 30) |
| `-l` | Convergence threshold for iterative parameter estimation (default: 0.001) |
| `-m` | Variance correction mode: `permute` (default) or `sample` |
| `-r` | Marker genotype file; required when `-m sample` is used |

#### Output

One file `<site>.middle` per splice site. The first line is tab-delimited and contains:

| Field | Description |
| ----- | ----------- |
| 1 | Model-fitting status: `Converged` or `Fixed` |
| 2 | Splice site name |
| 3 | Square root of the variance correction factor, √*r* |
| 4 | τ<sub>g</sub>, variance component of the GRM random effect (0 for `Fixed`) |

The following lines contain metacell-level quantities of the fitted null model (residuals, fitted usage ratios, total read counts and splice-supporting read counts), used in QTL mapping.

Collect the sites for which null models were built:

```bash
ls model_construct_QTL_mapping/model | sed -n 's/\.middle$//p' > model_construct_QTL_mapping/test_site.list
```

The metacell list for QTL mapping (`-m` below) must contain the metacells in the order used during model construction, i.e. the metacells of the phenotype file that are also present in the covariate file, in phenotype-file order, one per line.

### 3b. cis-sQTL mapping

Tests the association between each splice site and all variants within a cis window, using the pre-fitted null models. Process one chromosome at a time to enable parallelisation.

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
| `-o` | Output directory for QTL result files |
| `-c` | Chromosome to map (as named in the site list and genotype file) |
| `-v` | Genotype file in VCF/BCF format (must be indexed) |
| `-x` | Covariate (PC) file; must be the same file used for model construction |
| `-p` | Directory containing the pre-fitted null model files |
| `-w` | cis window half-width in bp (e.g., 500000 = ±500 kb around each splice site) |
| `-m` | Metacell list in the same order as the null model files (see above) |
| `-t` | P-value threshold; associations with P > threshold are not written |

#### Output

One file `<site>.result` per splice site, tab-delimited without header, with one line per tested variant:

| Column | Description |
| ------ | ----------- |
| 1 | Splice site |
| 2 | Variant ID |
| 3 | P value |
| 4 | Effect size (log-odds scale) |
| 5 | Standard error |

For sites with status `Fixed`, effect sizes are estimated under the fixed-effects model and therefore represent population-averaged (marginal) effects, which are typically attenuated relative to the conditional effects of the mixed model; P values are calibrated in both cases.

> **Note:** Please ensure all dependencies are installed before running the test. Refer to [Installation & Test Data](#installation--test-data) for details.

---

## Contact

- **Lab**: [Boxiang Liu Lab](https://github.com/boxiangliulab)
- **Issues**: Please report bugs via the [GitHub Issues](https://github.com/boxiangliulab/ISSAC/issues) page