---
name: issac
description: >
  Use this skill whenever the user asks about ISSAC (Integrative Single-cell Splicing Analysis
  and QTL Caller), including installation, running the pipeline, sQTL mapping, intron retention
  analysis, differential splicing, trans-sQTL, phenotype preparation, metacell detection, or
  interpreting ISSAC output files. Trigger whenever the user mentions ISSAC, single-cell sQTL,
  metacell splicing, SSUR, splice site usage ratio, scRNA-seq QTL mapping, GLMM splicing, or
  any of the ISSAC commands (junctools, juncstat, IR extract, pheno_group, pheno_output, model,
  QTL, DS, trans_QTL, IR_combine). Also trigger for questions about competitive introns vs
  intron retention in the context of single-cell data.
---

# ISSAC: Integrative Single-cell Splicing Analysis and QTL Caller

## What is ISSAC?

ISSAC maps **cis-sQTLs** (splicing quantitative trait loci) from single-cell RNA-seq data by:
- Aggregating cells into **metacells** to reduce sparsity
- Quantifying **splice site usage ratios (SSUR)** from junction and non-split reads
- Fitting a **generalized binomial linear mixed model (GLMM)** with GRM-based population structure correction

Two splicing phenotype types:
- **Competitive introns** — multiple competing splice sites within a cluster (classic sQTL)
- **Intron retention (IR)** — single intron clusters quantified via non-split reads

---

## Dependencies

```bash
conda install -c conda-forge gsl eigen nlopt openssl
conda install -c bioconda htslib=1.3
```

C++ libraries required: htslib 1.3, gsl, eigen3, nlopt, crypto (openssl)

---

## Installation

```bash
git clone --branch master https://github.com/boxiangliulab/ISSAC.git
cd ISSAC/build
rm -rf *
cmake ..
make
./ISSAC -h

# Set binary path for use throughout pipeline
ISSAC=path/to/ISSAC/build/ISSAC
```

---

## Pipeline Overview

```
BAM file (per metacell)
        │
        ▼
Step 1: Junction extraction + non-split read extraction
        │
        ▼
Step 2: Phenotype preparation (competitive introns & IR) + filtering
        │
        ▼
Step 3: Null GLMM fitting → cis-sQTL mapping
```

---

## Metacell Detection

Uses k-NN graph → Louvain clustering → merge small clusters (minimum size `m`).

Script: `analysis_pipeline/DLPFC_analysis/sQTL_mapping/metacell_calling.py`

Input: SCTransform-based PC embeddings from snRNA-seq gene expression matrix.

---

## Step 1: Junction & Non-split Read Extraction

### 1a. Junction Extraction

```bash
$ISSAC junctools extract \
  -a 8 \       # min reads supporting junction
  -m 50 \      # min intron length (bp)
  -M 500000 \  # max intron length (bp)
  -s RF \      # strand: RF for 10X 3' libraries; FR for 5' libraries
  $bamfile \
  -o $output_junc
```

Output `.junc` columns: chrom, read start, read end, CB-UMI-dedup count, strand, junction start, junction end, anchor lengths, cell barcode, UMI.

### 1b. Per-Metacell Junction Statistics

Aggregates junction counts per barcode using CB-UMI-based counting (removes PCR duplicates).

```bash
$ISSAC juncstat \
  -b $barcode \    # barcode list for this metacell
  -j $output_junc \
  -o $output_stat
```

Output `.stat`: `chrom:strand:start:end TAB count`

### 1c. Non-split Read Extraction (for IR)

```bash
$ISSAC IR extract \
  -s RF \          # strand specificity
  -b $barcode \    # metacell barcode list
  -t $site_list \  # intronic sites to quantify
  -a $bamfile \
  -o $output_nonsplit
```

---

## Step 2: Phenotype Preparation

### 2a. Competitive Intron Phenotype

```bash
$ISSAC pheno_group \
  -s $sample \           # file listing sample (metacell) names
  -j stat_file_dir/ \    # directory of .stat files
  -o output_prefix \
  -t 50 \                # min total reads per junction across all metacells
  -l log.out \
  -n 50 \                # min intron length (must match step 1)
  -x 500000              # max intron length (must match step 1)
```

Produces 4 output files:
- `.inclu_exclu` — included vs. excluded reads per splice site
- `.intron.out` — CB-UMI junction counts across samples
- `.refined` — per-intron counts across samples
- `.site` — per-site phenotype matrix `included:total` per sample

### 2b. Intron Retention Phenotype

```bash
$ISSAC IR_combine \
  -s $sample \
  -f nonsplit_file_dir/ \    # directory of .nonsplit files
  -l $single_intron_site \   # list of intronic sites to include
  -i $total_intron \         # .intron.out from step 2a
  -o output_prefix
```

### 2c. Phenotype Filtering

Apply to both competitive and IR phenotype files:

```bash
$ISSAC pheno_output \
  -r input.site \
  -o output.filtered \
  -p output.prop \
  -s 0.1 \    # min variance threshold
  -n 0.5      # max sparsity threshold (fraction missing)
```

`.filtered` — sites passing filters, format `included:total` per sample  
`.prop` — splice usage ratios (used for splice PC computation)

---

## Step 3: Model Construction & cis-sQTL Mapping

### GRM Preparation

```bash
# PLINK
plink --bfile pruned_output --make-grm-bin --out grm_output

# R: convert to text
library(plinkFile)
dat <- readGRM("path/to/grm")
write.table(as.data.frame(dat), "GRM.txt",
            sep = " ", row.names = TRUE, col.names = TRUE, quote = FALSE)
```

**Note:** Sample names in phenotype and PC files must follow `${donor}:${metacell_index}` format. `${donor}` must exist as row/column label in the GRM.

### 3a. Null Model Construction

```bash
$ISSAC model \
  -s $site_pheno \     # .filtered phenotype file
  -p $PC_file \        # genotype PC file
  -n 617 \             # number of individuals in GRM
  -g GRM.txt \
  -u output_model_dir/ \
  -t 10                # normalization iterations (×10)
```

After fitting, collect successfully built model sites:
```bash
ls model_dir/* | cut -d'/' -f3 | cut -d'.' -f1 > test_site.list
```

Model output rows: site name + params, residuals, π estimates, total counts, included counts.

### 3b. cis-sQTL Mapping

Run per-chromosome (parallelise across chromosomes):

```bash
$ISSAC QTL \
  -s $site_list \          # sites to test
  -o output_prefix \
  -c chr10 \               # one chromosome per job
  -v genotype.bcf \        # BCF format, must have .csi index
  -x $PC_file \
  -p model_dir/ \          # pre-fitted null model directory
  -w 500000 \              # cis window ±bp
  -m $common_name \        # sample name list matching model files
  -t 1                     # p-value threshold (discard if p > threshold)
```

Output `.result` columns: splice site, SNP (chrom:pos:ref:alt), p-value, effect size, SE.

---

## Differential Splicing

Reuses null model files directly (no re-fitting needed).

```bash
$ISSAC DS \
  -s site.list \
  -p model_dir/ \
  -m *.common \     # sample names
  -x *.PC \
  -g *.group \      # group labels: 0 = group A, 2 = group B (one per line)
  -o output_prefix
```

---

## Trans-sQTL Identification

```bash
$ISSAC trans_QTL \
  -s site.list \
  -p model_dir/ \
  -m *.common \
  -x *.PC \
  -c chr${i} \          # chromosome of splice sites
  -v chr${i}.recode.bcf \ # must have .csi index
  -w *.variant \        # list of variant IDs to test
  -o output_prefix \
  -t 1                  # p-value threshold
```

---

## Key Concepts

**SSUR (Splice Site Usage Ratio):** `included_reads / (included_reads + competing_or_nonsplit_reads)`. For IR events:
- `SSUR_a = ab / (ab + N_a)` at the donor site
- `SSUR_b = ab / (ab + N_b)` at the acceptor site

Where `ab` = junction reads spanning the intron, `N_a`/`N_b` = non-split reads at donor/acceptor.

**Imputation:** Sites with missing SSUR across samples are imputed using the mean of observed samples. A sparsity threshold of ≤50% missing is recommended (comparable to LeafCutter's ~40% default for pseudobulk).

**Working weights:** The GLMM uses binomial working weights. When N is imputed (unreliable), weighting by N·p(1-p) can introduce spurious precision for imputed observations; using p(1-p) alone is a valid conservative alternative when N is homogeneous or unreliable.

**Strand flags:**
- `RF` → 10X 3' libraries (most common)
- `FR` → 10X 5' libraries
