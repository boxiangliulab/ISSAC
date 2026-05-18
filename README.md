# ISSAC: Integrative Single-cell Splicing Analysis and QTL Caller

### ISSAC is a scalable tool for single-cell sQTL mapping by modeling metacell splice site usage ratios with a generalized binomial mixed model.

![image](https://github.com/tibettiger/ISSAC/blob/ISSAC_analysis/img/figure1_v2.png)

---

## Table of Contents

- [Dependencies](#dependencies)
- [Installation](#installation)
- [Pipeline Overview](#pipeline-overview)
- [Metacell Detection](#metacell-detection)
- [Step 1: Junction & Non-split Read Extraction](#step-1-junction--non-split-read-extraction)
- [Step 2: Phenotype Preparation](#step-2-phenotype-preparation)
- [Step 3: Model Construction & cis-sQTL Mapping](#step-3-model-construction--cis-sqtl-mapping)
- [Differential Splicing](#differential-splicing)
- [Trans-sQTL Identification](#trans-sqtl-identification)

---

## Dependencies

The following C++ libraries must be available in your environment before building ISSAC:

1. htslib 1.3
2. gsl
3. eigen3
4. nlopt
5. crypto (openssl)

```bash
conda install -c conda-forge gsl eigen nlopt openssl
conda install -c bioconda htslib=1.3
```

---

## Installation

```bash
git clone --branch master https://github.com/boxiangliulab/ISSAC.git
cd ISSAC/build
rm -rf *
cmake ..
make
./ISSAC -h
```

A successful installation will print:

```
Usage:          ISSAC <command> [options]
Command:        Integrative single-cell splicing analysis and QTL caller
```

Set the path to the compiled binary for use throughout the pipeline:

```bash
ISSAC=path/to/ISSAC/build/ISSAC
```

---

## Pipeline Overview

```
BAM file
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

To improve statistical power while keeping sparsity low, ISSAC operates on metacells rather than individual cells. For each donor:

1. PC embeddings from the snRNA-seq gene expression matrix (SCTransformed-based) are used to build a k-nearest-neighbours graph.
2. Louvain clustering is applied with `m` as the minimum cluster size.
3. Clusters smaller than `m` cells are iteratively merged into adjacent larger clusters.
4. Each final cluster of ≥ `m` cells is treated as one metacell for sQTL mapping.

The metacell detection script is available at:
[metacell_calling.py](https://github.com/boxiangliulab/ISSAC/blob/ISSAC_analysis/analysis_pipeline/DLPFC_analysis/sQTL_mapping/metacell_calling.py)

---

## Step 1: Junction & Non-split Read Extraction

Extract splicing junctions and non-split reads from BAM files. These provide raw per-barcode evidence for splicing events and intron retention. Repeat for each metacell BAM file; a corresponding `.bai` index file must be present alongside each BAM.

### 1a. Junction Extraction

Extracts splice junctions from a 10X scRNA-seq BAM file.

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
| `-m` | Minimum intron length (bp); junctions spanning smaller distances are excluded |
| `-M` | Maximum intron length (bp); junctions spanning larger distances are excluded |
| `-s` | Strand specificity (`RF` for 10X 3′ libraries; `FR` for 5′ libraries) |
| `-o` | Output juncfiles' prefix |

**Example output (`.junc`):**

```
chr1    960660  961744  6       +       960800  961628  140,116 TGCACCTAGGTCGGAT        GCGCGCGCGT
chr1    961651  961879  4       +       961750  961825  99,54   TGCACCTAGGTCGGAT        GCGCGCGCGT
chr1    1013487 1014125 17      +       1013576 1013983 89,142  TGGGCGTAGTACGATA        CTTCACAGAT
```

Columns: (1) chrom, (2) read start, (3) read end, (4) CB-UMI–deduplicated junction read count, (5) strand, (6) junction start, (7) junction end, (8) anchor lengths, (9) cell barcode, (10) UMI.

### 1b. Per-Metacell Junction Statistics

Aggregates junction read counts per barcode (metacell) using CB-UMI–based counting to avoid PCR amplification bias.

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
| `-o` | Output statistics file (CB-UMI–based counts per junction per metacell) |

**Example output (`.stat`):**

```
chr10:+:124484326:124488421     1
chr10:+:124801901:124806830     1
chr10:+:125719915:125721203     3
chr10:+:130136512:130145210     1
chr10:+:132397351:132404625     2
```

Columns: (1) intron junction (chrom:strand:start:end), (2) CB-UMI–based read count.

### 1c. Non-split Read Extraction

Extracts reads that do **not** span a splice junction, used for intron retention (IR) quantification at sites of interest.

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
| `-t` | List of intronic sites to quantify non-split read coverage |
| `-a` | Input BAM file |
| `-o` | Output non-split read file |

**Example output (`.nonsplit`):**

```
chr10:100523886:-:CATCCCACATTGTAGC:CGCGGTGGCGGT 5
chr10:100523886:-:GTAGCTAAGACTCTTG:ATTGCCTTAATT 1
chr10:100523886:-:TCCTCTTTCATCGGGC:ACAGGTAGATCT 1
chr10:100525399:-:TGGCGTGTCAGCGCAC:GGACTTAGCAAG 1
```

Each line: non-split CB-UMI–based read counts crossing a splice site (chrom:pos:strand:barcode:UMI).

---

## Step 2: Phenotype Preparation

Constructs per-metacell splicing phenotype matrices from junction and non-split read files. Two types of splicing events are handled:

- **(a) Competitive introns** — splice sites with detectable competing introns (classic sQTL signal)
- **(b) Single intron clusters** — intron retention events quantified from non-split reads

### 2a. Competitive Intron Phenotype Preparation

Groups splice sites across metacells and generates per-site usage ratios

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
| `-t` | Minimum total reads per junction across all metacells |
| `-l` | Log file for intermediate results |
| `-n` | Minimum intron length (bp); must match junction extraction |
| `-x` | Maximum intron length (bp); must match junction extraction |

`pheno_group` produces four output files:

**`.inclu_exclu`** — for each splice site: the supporting (included) and competing (excluded) intron reads.

```
chr10:+:73 110985627 included 110919657:110985627 excluded 110919657:110951607 110919657:110964124
chr10:+:76 111114416 included 111114416:111114902 111114416:111117959 excluded
```

**`.intron.out`** — CB-UMI–based junction read counts for each intron across all samples, the same sequence as in the sample list.

```
chr10:+:233 chr10:+:22925636:22931995 2 0 2 0 1 0
chr10:+:234 chr10:+:22932044:22946143 1 1 1 0 0 0
chr10:+:234 chr10:+:22946261:22955806 1 0 0 0 0 0
```

**`.refined`** — each intron's CB-UMI–based junction reads across all samples.

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

### 2b. Single Intron Cluster (Intron Retention) Phenotype Preparation

Combines per-metacell non-split read files into a site-level phenotype matrix, quantifying intron retention as the ratio of split reads to total reads (including both split and non-split reads) at each site.

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
| `-f` | Directory containing per-metacell `.nonsplit` files |
| `-l` | List of intronic sites to include |
| `-i` | File for total intron read counts across all metacells (`.intron.out` from 2a) |
| `-o` | Output prefix for IR phenotype files |

### 2c. Phenotype Filtering

Filters out splice sites with low variance or high sparsity, retaining only informative sites for QTL mapping. Apply to both competitive intron and IR phenotype files.

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
| `-p` | Output file for per-site usage proportions (`.prop`) |
| `-s` | Minimum variance threshold; sites below this are excluded |
| `-n` | Maximum sparsity threshold; sites exceeding this fraction of missing values are excluded |

**Example `.filtered` output** (header row lists samples; subsequent rows are sites that passed filtering):

```
JP_RIK_H002:1 JP_RIK_H017:1 JP_RIK_H019:1 JP_RIK_H024:1 JP_RIK_H026:1
chr10:+:104254499 1:2 2:13 2:15 0:6 1:7
chr10:+:104254573 0:1 0:2 0:2 1:2 1:2
chr10:+:104255162 2:2 13:13 15:15 6:6 7:7
```

**Example `.prop` output** (splice site usage ratios; used for splice PC computation):

```
chr10:-:120793306 0.000000 0.000000 0.000000 0.000000 0.000000
chr10:-:129036456 0.666667 0.666667 1.000000 1.000000 0.666667
chr10:-:16816972  0.846154 1.000000 0.777778 0.888889 1.000000
```

---

## Step 3: Model Construction & cis-sQTL Mapping

Fits a binomial mixed model (GLMM) per splice site using a genetic relatedness matrix (GRM) to control for population structure, then performs cis-sQTL mapping within a defined window.

### GRM Preparation

Generate a GRM from genotype data (pruned genotype) using PLINK or GCTA:

```bash
plink --bfile pruned_output --make-grm-bin --out grm_output
```

Convert the binary GRM to a text file for ISSAC using R:

```r
library(plinkFile)

dat <- readGRM("path/to/grm")
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

### PC File Format

Note: all sample names in the phenotype and PC files must follow the format `${donor}:${metacell_index}`, and `${donor}` must exist as a row/column label in the GRM file.

**Example `.PC` file:**

```
JP_RIK_H002:1   JP_RIK_H002:3   JP_RIK_H003:1   JP_RIK_H003:10
-2.78345826143274   -2.54253101018537   -1.28028808817251   -1.65413653429238
1.36534270723467    0.625185320131184   -0.747500429552883   0.229720830781539
```

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
| `-s` | Filtered splicing phenotype file (`.filtered`) |
| `-p` | Genotype PC file for covariate correction |
| `-n` | Number of individuals in the GRM file |
| `-g` | GRM file (`.txt`) |
| `-u` | Output directory/prefix for fitted null model files |
| `-t` | Number of normalization parameter estimation iterations (×10) |

**Example model file output:**

```
Splice_site     chr19:+:34254614        0.672887        2.04307e-14
residuals       0.031697        0.0262904       0.015532        0.030618
pi              0.984151        0.986855        0.984468        0.984691
total           2       2       1       2
y               2       2       1       2
```

Row meanings: (1) site name, normalization parameter, variance component of random effect; (2) per-sample null residuals; (3) per-sample null π estimates; (4) total CB-UMI counts per sample; (5) CB-UMI counts supporting site usage per sample. Sample order matches the PC and phenotype files.

Collect sites for which null models were successfully built:

```bash
ls model_construct_QTL_mapping/model/* \
  | cut -d '/' -f 3 \
  | cut -d '.' -f 1 \
  > model_construct_QTL_mapping/test_site.list
```

### 3b. cis-sQTL Mapping

Tests association between each splice site and all SNPs within a cis window using pre-fitted null models. Process one chromosome at a time to enable parallelisation.

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
| `-x` | PC file (must match covariates used during model construction) |
| `-p` | Directory containing pre-fitted null model files |
| `-w` | cis window size in bp (SNPs within ±`w` bp of each site are tested) |
| `-m` | File listing sample names in the same order as the null model files |
| `-t` | P-value output threshold; associations with p > threshold are discarded |

**Example `.result` output:**

```
chr6:-:100847625        chr6:100351323:C:T      0.00253261      0.0413224       0.0136855
chr6:-:100847625        chr6:100352062:G:A      0.882403        0.00847691      0.0573062
chr6:-:100847625        chr6:100352246:G:A      0.451096        -0.00176255     0.00233887
```

Columns: (1) splice site, (2) SNP (chrom:pos:ref:alt), (3) p-value, (4) effect size, (5) standard error of effect size.

---

## Differential Splicing

Null model files and metacell group labels can be used directly for differential splicing analysis without re-fitting models.

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
| `-s` | Site list (from model construction step) |
| `-p` | Directory containing pre-fitted null model files |
| `-m` | File listing sample names matching the phenotype and PC files |
| `-x` | PC file |
| `-g` | Group label file (one label per line: `0` for group A, `2` for group B) |
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

Tests associations between splice site usage ratios and SNPs located distally from the splice site.

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
| `-s` | Site list (from model construction step) |
| `-p` | Directory containing pre-fitted null model files |
| `-m` | File listing sample names |
| `-x` | PC file |
| `-c` | Chromosome of the splice sites |
| `-v` | Genotype BCF file (must have a `.csi` index) |
| `-w` | File listing variant IDs to test |
| `-o` | Output prefix |
| `-t` | P-value output threshold; associations with p > threshold are discarded |




