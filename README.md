# ISSAC: Integrative Single-cell Splicing Analysis and QTL Caller

### ISSAC is a scalable tool for single-cell sQTL mapping that models metacell-level splice site usage with a binomial generalized linear mixed model (GLMM).

![ISSAC overview](https://github.com/tibettiger/ISSAC/blob/ISSAC_analysis/img/figure1_v4.png)

---

## Table of Contents

- [Installation & Test Data](#installation--test-data)
  - [Option 1: Native Linux build](#option-1-native-linux-build)
  - [Option 2: Docker (Mac, Windows, Linux)](#option-2-docker-mac-windows-linux)
  - [Option 3: Apptainer/Singularity (HPC users)](#option-3-apptainersingularity-hpc-users)
- [Notes](#notes)
- [Pipeline Overview](#pipeline-overview)
- [Metacell Detection](#metacell-detection)
- [Step 1: Junction & Non-split Read Extraction](#step-1-junction--non-split-read-extraction)
- [Step 2: Phenotype Preparation](#step-2-phenotype-preparation)
- [Step 3: Model Construction & cis-sQTL Mapping](#step-3-model-construction--cis-sqtl-mapping)
- [Differential Splicing](#differential-splicing)
- [Trans-sQTL Identification](#trans-sqtl-identification)

---

## Installation & Test Data

ISSAC can be run in three ways: natively on Linux (after building from source), via Docker (Mac, Windows and Linux), or via Apptainer/Singularity (HPC clusters). A small example dataset is provided under `test_data/` to verify your installation and explore ISSAC's functionality.

### Clone the repository

```bash
git clone -b master https://github.com/boxiangliulab/ISSAC.git
cd ISSAC
```

### Option 1: Native Linux build

#### Dependencies

The following libraries must be available before compiling ISSAC:

1. htslib
2. gsl
3. eigen3
4. nlopt
5. crypto (openssl)

Build tools: CMake (≥ 3.10), make and a C++11 compiler.

The simplest way is to install everything into a dedicated conda environment:

```bash
conda create -n issac_env -y -c conda-forge -c bioconda \
    gsl eigen nlopt openssl htslib zlib cmake make compilers
conda activate issac_env
```

#### Build

From the root of the repository:

```bash
mkdir -p build && cd build
cmake .. -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
make -j4
cd ..
```

The executable is written to `build/ISSAC`. The commands in this README call ISSAC through the variable `$ISSAC`:

```bash
export ISSAC=$(pwd)/build/ISSAC
$ISSAC -h
```

#### Run the test pipeline

```bash
cd test_data/
bash linux_test.sh
```

### Option 2: Docker (Mac, Windows, Linux)

No local build is required. This pulls a prebuilt image supporting both `amd64` and `arm64` (Apple Silicon) architectures.

```bash
cd test_data/
bash docker_test.sh
```

> **Windows users:** run this script inside Git Bash or WSL, not PowerShell/CMD.

A successful run shows intermediate output such as `Converged` and `Genotype read in successfully`, and the pipeline proceeds through `Start pvalue computing!` without errors.

To run ISSAC manually rather than via the wrapper script:

```bash
docker pull yuntian1999/issac:v1.3
docker run --rm yuntian1999/issac:v1.3 -h
```

To use your own data, mount your working directory and set it as the container's working directory (the image's entrypoint is `ISSAC`, so pass the subcommand directly):

```bash
docker run --rm -v "$PWD":"$PWD" -w "$PWD" yuntian1999/issac:v1.3 QTL [options]
```

### Option 3: Apptainer/Singularity (HPC users)

```bash
module load apptainer
cd test_data/
bash apptainer_test.sh
```

To run ISSAC manually:

```bash
apptainer pull issac.sif docker://yuntian1999/issac:v1.3
apptainer exec issac.sif ISSAC -h
```

To use your own data, bind mount your working directory:

```bash
apptainer exec --bind "$PWD":"$PWD" --pwd "$PWD" issac.sif ISSAC QTL [options]
```

When using a container, replace `$ISSAC` in the commands below with `docker run --rm -v "$PWD":"$PWD" -w "$PWD" yuntian1999/issac:v1.3` or `apptainer exec --bind "$PWD":"$PWD" --pwd "$PWD" issac.sif ISSAC`.

---

## Notes

- Docker images are published at [`yuntian1999/issac`](https://hub.docker.com/r/yuntian1999/issac) and support both `linux/amd64` and `linux/arm64`.
- The Apptainer `.sif` container is built directly from the same Docker image, so the Docker and Apptainer installations run identical binaries.

---

## Pipeline Overview

```
snRNA-seq data ──► Metacell detection (per donor)
        │
        ▼
BAM file + metacell barcode lists
        │
        ▼
Step 1: Junction extraction + non-split read extraction
        │
        ▼
Step 2: Phenotype preparation (competitive introns & intron retention) + filtering
        │
        ▼
Step 3: Null GLMM fitting → cis-sQTL mapping
```

---

## Metacell Detection

To improve statistical power and reduce sparsity, ISSAC operates on metacells rather than individual cells. Metacells are constructed separately for each donor, so that every metacell contains cells from a single donor:

1. PC embeddings of the gene expression matrix (log-normalized or SCTransform-normalized) are used to build a k-nearest-neighbour graph (k = `m`).
2. Louvain clustering is applied to the graph.
3. Clusters with fewer than `m` cells are dissolved, and their cells are reassigned to the nearest retained cluster (by Euclidean distance to the cluster centroid in PC space).
4. Each final cluster (≥ `m` cells) is treated as one metacell. Donors with fewer than 2`m` cells, or without at least two clusters of ≥ `m` cells, are represented by a single metacell.

The metacell detection script is available at:
[metacell_calling.py](https://github.com/boxiangliulab/ISSAC/blob/ISSAC_analysis/analysis_pipeline/DLPFC_analysis/sQTL_mapping/metacell_calling.py)

Metacell names follow the format `${donor}:${metacell_index}`, which is required by the downstream steps.

---

## Step 1: Junction & Non-split Read Extraction

Extract splice junctions and non-split reads from BAM files. These provide raw per-barcode evidence for splicing events and intron retention. Junction extraction (1a) is run once per BAM file; junction statistics (1b) and non-split read extraction (1c) are run once per metacell, using the barcode list of that metacell. A `.bai` index file must be present alongside each BAM.

### 1a. Junction Extraction

Extracts splice junctions from a 10x Genomics scRNA-seq BAM file.

```bash
bamfile=junctions_nonsplit_extract/test.bam
output_junc=junctions_nonsplit_extract/test.junc

$ISSAC junctools extract \
  -a 8 \
  -m 50 \
  -M 500000 \
  -s RF \
  -o $output_junc \
  $bamfile
```

| Flag | Description |
|------|-------------|
| `-a` | Minimum anchor length (bp) |
| `-m` | Minimum intron length (bp); shorter junctions are excluded |
| `-M` | Maximum intron length (bp); longer junctions are excluded |
| `-s` | Strand specificity (`RF` for 10x 3′ libraries; `FR` for 5′ libraries) |
| `-o` | Output junction file (`.junc`) |

**Example output (`.junc`):**

```
chr1    960660  961744  6       +       960800  961628  140,116 TGCACCTAGGTCGGAT        GCGCGCGCGT
chr1    961651  961879  4       +       961750  961825  99,54   TGCACCTAGGTCGGAT        GCGCGCGCGT
chr1    1013487 1014125 17      +       1013576 1013983 89,142  TGGGCGTAGTACGATA        CTTCACAGAT
```

Columns: (1) chromosome, (2) read start, (3) read end, (4) number of junction-spanning reads sharing the same CB–UMI, (5) strand, (6) junction (intron) start, (7) junction (intron) end, (8) anchor lengths, (9) cell barcode (CB), (10) UMI.

### 1b. Per-Metacell Junction Statistics

Aggregates junction reads over all barcodes belonging to one metacell. Reads are counted by unique CB–UMI to avoid PCR amplification bias.

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
| `-b` | Barcode list file defining the cells that belong to this metacell |
| `-j` | Junction file produced in 1a |
| `-o` | Output statistics file (CB–UMI-based counts per junction for this metacell) |

**Example output (`.stat`):**

```
chr10:+:124484326:124488421     1
chr10:+:124801901:124806830     1
chr10:+:125719915:125721203     3
chr10:+:130136512:130145210     1
chr10:+:132397351:132404625     2
```

Columns: (1) intron (chrom:strand:start:end), (2) CB–UMI-based read count.

### 1c. Non-split Read Extraction

Extracts reads that do **not** span a splice junction but cover the splice sites of interest; these are used for intron retention (IR) quantification.

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
| `-s` | Strand specificity (`RF` for 3′; `FR` for 5′) |
| `-b` | Barcode list for the metacell |
| `-t` | List of splice sites at which non-split read coverage is quantified |
| `-a` | Input BAM file |
| `-o` | Output non-split read file |

**Example output (`.nonsplit`):**

```
chr10:100523886:-:CATCCCACATTGTAGC:CGCGGTGGCGGT 5
chr10:100523886:-:GTAGCTAAGACTCTTG:ATTGCCTTAATT 1
chr10:100523886:-:TCCTCTTTCATCGGGC:ACAGGTAGATCT 1
chr10:100525399:-:TGGCGTGTCAGCGCAC:GGACTTAGCAAG 1
```

Each line: a splice site and CB–UMI (chrom:pos:strand:barcode:UMI), followed by the number of non-split reads covering that site.

---

## Step 2: Phenotype Preparation

Constructs metacell-level splicing phenotype matrices from the junction and non-split read files. Two types of splicing events are handled:

- **(a) Competitive introns** — splice sites used by two or more competing introns; usage is quantified from split reads.
- **(b) Single-intron sites** — splice sites used by a single intron; usage is quantified against non-split reads (intron retention).

### 2a. Competitive Intron Phenotype Preparation

Groups introns sharing splice sites across metacells and generates per-site usage counts.

```bash
ls splice_phenotype_prepare/stat_file/ \
  | grep '\.stat$' \
  | sed 's/\.stat$//' \
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
| `-t` | Minimum total reads per intron across all metacells |
| `-l` | Log file for intermediate results |
| `-n` | Minimum intron length (bp); should match 1a |
| `-x` | Maximum intron length (bp); should match 1a |

`pheno_group` produces four output files:

**`.inclu_exclu`** — for each splice site: the introns using the site (included) and the competing introns (excluded).

```
chr10:+:73 110985627 included 110919657:110985627 excluded 110919657:110951607 110919657:110964124
chr10:+:76 111114416 included 111114416:111114902 111114416:111117959 excluded
```

**`.intron.out`** — CB–UMI-based junction read counts for each intron across all samples, in the same order as the sample list.

```
chr10:+:233 chr10:+:22925636:22931995 2 0 2 0 1 0
chr10:+:234 chr10:+:22932044:22946143 1 1 1 0 0 0
chr10:+:234 chr10:+:22946261:22955806 1 0 0 0 0 0
```

**`.refined`** — introns in each intron cluster, with their total CB–UMI-based junction reads across all samples (start:end:reads).

```
chr10:+:111 11852215:11862911:125 11852215:11866530:640
chr10:+:112 118730410:118754395:128
chr10:+:115 119207969:119326515:530 119207969:119380814:38 119326611:119380814:182
```

**`.site`** — per-site phenotype matrix: `included_reads:total_reads` for each sample.

```
chr10:-:16782084 13:13 11:11 17:18 9:9 6:6 15:15 15:15
chr10:-:16816972 11:13 11:11 14:18 8:9 6:6 15:15 14:15
chr10:-:16817084 9:11 12:12 13:17 8:9 7:7 13:13 15:16
```

### 2b. Single-Intron Site (Intron Retention) Phenotype Preparation

Combines the per-metacell non-split read files into a site-level phenotype matrix, quantifying usage at each site as split reads / (split reads + non-split reads).

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
| `-s` | Sample list file (same as 2a) |
| `-f` | Directory containing per-metacell `.nonsplit` files |
| `-l` | List of single-intron sites to include |
| `-i` | Intron read counts across all metacells (`.intron.out` from 2a) |
| `-o` | Output prefix for IR phenotype files |

### 2c. Phenotype Filtering

Removes splice sites with low variability or high sparsity, retaining informative sites for QTL mapping. Apply to both the competitive intron and the IR phenotype files.

**Competitive intron sites:**

```bash
s=0.1   # minimum variance threshold
n=0.5   # maximum sparsity threshold

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
| `-r` | Input site phenotype file (`.site`) |
| `-o` | Filtered output phenotype file (`.filtered`) |
| `-p` | Output file of per-site usage proportions (`.prop`) |
| `-s` | Minimum variance threshold; sites below this are excluded |
| `-n` | Maximum sparsity threshold; sites with a larger fraction of missing values are excluded |

**Example `.filtered` output** (the header row lists samples; subsequent rows are sites that passed filtering):

```
JP_RIK_H002:1 JP_RIK_H017:1 JP_RIK_H019:1 JP_RIK_H024:1 JP_RIK_H026:1
chr10:+:104254499 1:2 2:13 2:15 0:6 1:7
chr10:+:104254573 0:1 0:2 0:2 1:2 1:2
chr10:+:104255162 2:2 13:13 15:15 6:6 7:7
```

**Example `.prop` output** (splice site usage ratios; used to compute splicing PCs):

```
chr10:-:120793306 0.000000 0.000000 0.000000 0.000000 0.000000
chr10:-:129036456 0.666667 0.666667 1.000000 1.000000 0.666667
chr10:-:16816972  0.846154 1.000000 0.777778 0.888889 1.000000
```

---

## Step 3: Model Construction & cis-sQTL Mapping

For each splice site, ISSAC fits a binomial GLMM with two random effects: a donor-level genetic random effect, whose covariance is given by the metacell-level genetic relationship matrix (GRM; this accounts for the correlation among metacells from the same donor and for relatedness between donors), and an observation-level random effect that captures overdispersion. Population structure and other confounders are adjusted for through covariates (e.g., genotype PCs and splicing PCs). cis-sQTL mapping is then performed with a score test within a defined window.

### GRM Preparation

Generate a GRM from LD-pruned genotype data using PLINK (v1.9) or GCTA:

```bash
plink --bfile pruned_output --make-grm-bin --out grm_output
# or: gcta64 --bfile pruned_output --make-grm --out grm_output
```

Convert the binary GRM to a text file for ISSAC using R:

```r
library(plinkFile)

dat <- readGRM("grm_output")   # prefix of the .grm.bin/.grm.id files
write.table(as.data.frame(dat), "GRM.txt",
            sep = " ", row.names = TRUE, col.names = TRUE, quote = FALSE)
```

**Example `GRM.txt`:**

```
JP_RIK_H001 JP_RIK_H002 JP_RIK_H003
JP_RIK_H001  0.989882    0.030982    0.030328
JP_RIK_H002  0.030982    0.982672    0.011874
JP_RIK_H003  0.030328    0.011874    0.986594
```

### Covariate (PC) File Format

The covariate file contains one column per metacell and one row per covariate (e.g., genotype PCs, splicing PCs computed from the `.prop` file, and other covariates such as sex or age). The first row lists the sample names.

All sample names in the phenotype and covariate files must follow the format `${donor}:${metacell_index}`, and each `${donor}` must appear as a row/column label in the GRM file.

**Example `.PC` file:**

```
JP_RIK_H002:1   JP_RIK_H002:3   JP_RIK_H003:1   JP_RIK_H003:10
-2.78345826143274   -2.54253101018537   -1.28028808817251   -1.65413653429238
1.36534270723467    0.625185320131184   -0.747500429552883   0.229720830781539
```

### 3a. Null Model Construction

Pre-fits the null GLMM (without genotype) for each splice site, so that it does not need to be refitted for every SNP during QTL mapping. In this step, ISSAC also estimates the site-specific variance ratio *r*, which corrects the discrepancy between the true variance of the score statistic and its variance under the working-weight matrix W. Two estimation modes are available:

| Mode | Option | Recommended for |
|------|--------|-----------------|
| **Permutation** (default) | `-m permute` | Cohorts with limited between-donor relatedness (most population-based single-cell cohorts) |
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

Because P depends on the estimated variance components, this mode requires a stably converged null model. Splice sites are handled according to their model-fitting status:

| Model-fitting status | Variance correction |
|----------------------|---------------------|
| Converged, τ<sub>g</sub> > 10<sup>-6</sup> | SAIGE-style G<sup>T</sup>PG / G<sup>T</sup>WG ratio |
| Not converged, or converged to the lower boundary (τ<sub>g</sub> ≤ 10<sup>-6</sup>) | Site excluded |

**Marker genotype file (`-r`).** Tab-delimited; the first column is the donor ID and each remaining column is one marker with genotype dosages (0/1/2; missing values as `NA`). An optional header line is allowed. **Donors must be in the same order as in the GRM file** (donor IDs are not matched by name).

```
donor_id   rs1001   rs1002   rs1003   ...
D001       0        1        2
D002       1        0        1
D003       2        1        0
```

We recommend at least 30 common (MAF ≥ 0.05), approximately independent markers located on chromosomes other than those of the tested splice sites, so that they have no *cis* effect on splicing.

#### Options

| Flag | Description |
|------|-------------|
| `-s` | Filtered splicing phenotype file (`.filtered`) |
| `-p` | Covariate file (genotype PCs, splicing PCs and other covariates; see format above) |
| `-n` | Number of individuals in the GRM file |
| `-g` | GRM file (`.txt`) |
| `-u` | Output directory for the fitted null model files |
| `-t` | Number of permutations for estimating *r*, in units of 100 (default: 10, i.e., 1,000 permutations). Used only with `-m permute` |
| `-v` | GRM sparsification threshold; relatedness values below this threshold are set to zero |
| `-i` | Maximum number of iterations for estimating fixed and random effects |
| `-l` | Convergence threshold for iterative parameter estimation (default: 0.001) |
| `-m` | Variance ratio estimation mode: `permute` (default) or `sample` |
| `-r` | Marker genotype file; required when `-m sample` is used |

**Example model file output:**

```
Converged     chr19:+:34254614        0.89        0.015
residuals       0.031697        0.0262904       0.015532        0.030618
pi              0.984151        0.986855        0.984468        0.984691
total           2       2       1       2
y               2       2       1       2
```

Row meanings:
1. Model-fitting status (`Converged`, or `Fixed` when the mixed model did not converge and ISSAC fell back to a fixed-effects model), site name, variance ratio *r*, and the estimated genetic variance component τ<sub>g</sub> (reported as −1 for `Fixed` sites; values ≤ 10<sup>-6</sup> indicate convergence to the lower boundary);
2. per-sample null residuals;
3. per-sample null π estimates;
4. total CB–UMI counts per sample;
5. CB–UMI counts supporting site usage per sample.

The sample order matches the PC and phenotype files.

Collect the sites for which null models were successfully built:

```bash
ls model_construct_QTL_mapping/model/ \
  | cut -d '.' -f 1 \
  > model_construct_QTL_mapping/test_site.list
```

### 3b. cis-sQTL Mapping

Tests the association between each splice site and all SNPs within its cis window using the pre-fitted null models. Run one chromosome per job to parallelise.

```bash
genotype=model_construct_QTL_mapping/test.bcf   # must have a .csi index
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
| `-c` | Chromosome to map (run one chromosome per job for parallelisation) |
| `-v` | Genotype file in BCF format (must have a `.csi` index) |
| `-x` | Covariate (PC) file; must be identical to the one used in 3a |
| `-p` | Directory containing the pre-fitted null model files |
| `-w` | cis window size in bp (SNPs within ±`w` bp of each site are tested) |
| `-m` | File listing sample (metacell) names in the same order as the phenotype and null model files |
| `-t` | P-value output threshold; associations with p > threshold are not written |

**Example `.result` output:**

```
chr6:-:100847625        chr6:100351323:C:T      0.00253261      0.0413224       0.0136855
chr6:-:100847625        chr6:100352062:G:A      0.882403        0.00847691      0.0573062
chr6:-:100847625        chr6:100352246:G:A      0.451096        -0.00176255     0.00233887
```

Columns: (1) splice site, (2) SNP (chrom:pos:ref:alt), (3) p-value, (4) effect size (β, on the logit scale), (5) standard error of the effect size.

> For sites whose null model converged to the lower boundary (τ<sub>g</sub> ≤ 10<sup>-6</sup>) or fell back to a fixed-effects model, the p-value is still calibrated by *r*, but the effect size is computed without the variance ratio and should be interpreted with caution. These sites can be identified from the first row of the corresponding null model file.

---

## Differential Splicing

The null model files and metacell group labels can be used directly for differential splicing analysis without refitting the models.

```bash
$ISSAC DS \
  -s site.list \
  -p $model_file_pos \
  -m *.common \
  -x *.PC \
  -g *.group \
  -o $output_pos
```

| Flag | Description |
|------|-------------|
| `-s` | Site list (from the model construction step) |
| `-p` | Directory containing the pre-fitted null model files |
| `-m` | File listing sample names, in the same order as the phenotype and PC files |
| `-x` | Covariate (PC) file |
| `-g` | Group label file: one label per line, in the same sample order as `-m` (`0` for group A, `2` for group B) |
| `-o` | Output prefix |

**Example `.group` file:**

```
0
0
2
2
0
2
```

---

## Trans-sQTL Identification

Tests associations between splice site usage and SNPs located distally from the splice site.

```bash
$ISSAC trans_QTL \
  -s site.list \
  -p $model_file_pos \
  -m *.common \
  -x *.PC \
  -c chr${i} \
  -v chr${i}.recode.bcf \
  -w *.variant \
  -o $output_pos \
  -t 1
```

| Flag | Description |
|------|-------------|
| `-s` | Site list (from the model construction step) |
| `-p` | Directory containing the pre-fitted null model files |
| `-m` | File listing sample names |
| `-x` | Covariate (PC) file |
| `-c` | Chromosome of the variants being tested |
| `-v` | Genotype BCF file for that chromosome (must have a `.csi` index) |
| `-w` | File listing variant IDs to test |
| `-o` | Output prefix |
| `-t` | P-value output threshold; associations with p > threshold are not written |





