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

## Test data
```bash
#
git clone https://github.com/boxiangliulab/ISSAC.git
cd test_data/
bash test_detail.sh

---


---


---

## License

---

## Contact
- **Lab**: [Boxiang Liu Lab](https://github.com/boxiangliulab)
- **Issues**: Please report bugs via the [GitHub Issues](https://github.com/boxiangliulab/ISSAC/issues) page
