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

## Components
| File | Description |
|---|---|
| `junctions_extract.cpp` | Single-cell junction extraction (adapted from regtools) |
| `junctions_extract.h` | Header file for junction extraction module |

> `junctions_extract.cpp` and `junctions_extract.h` were adapted from [regtools](https://github.com/griffithlab/regtools) to enable junction extraction at the single-cell level.

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

> **Note:** Please ensure all dependencies are installed before running the test. Refer to the [Installation](#installation) section for details.

---

## License

---

## Contact
- **Lab**: [Boxiang Liu Lab](https://github.com/boxiangliulab)
- **Issues**: Please report bugs via the [GitHub Issues](https://github.com/boxiangliulab/ISSAC/issues) page
