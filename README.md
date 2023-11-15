# ecDNAfinder

**Last updated: November 14, 2023**

## Prerequisites

### Packages

* **GNU Parallel** (tested on 20230722)
* **Python** (tested on 3.12.0)
  * **pygini** (tested on 1.0.1)
  * **pandas** (tested on 2.1.2)
  * **numpy** (tested on 1.26.0)
* **R** (tested on 4.3.2)
  * **CMplot** (tested on 4.4.3)
* (other packages required as dependencies)

### Data

* **Single-cell Hi-C contact matrices**
* **Copy number variation data**
* **Trained linear model**

Note: see **File format requirements** for examples. File names are currently hard-coded in `ecDNAfinder` script. Modify them if needed.

## Installation

We recommend to run `ecDNAfinder` in a conda/mamba environment.
```bash
conda create -n ecDNAfinder python parallel pandas numpy
pip install pygini
```
Then clone this repository to a local directory.

## Usage

### 1. ecDNAfinder

**Example:**

```bash
sh ecDNAfinder -i DATA -o OUT -p 0.95 -t 16

sh ecDNAfinder -i DATA -o OUT -p 0.95 -t 16 -s true
```

**Arguments:**

`-i` input directory. See **Directory hierarchy requirements** for details.

`-o` output directory. See **Directory hierarchy requirements** for details.

`-p` cutoff threshold. Value represents probability from the logistic regression.

`-t` number of threads. Current version occasionally crashes if the number of threads is too large.

`-s` (optional) summary-only mode. Use this option to skip redundant cell processing.

**Output example:**

```
> ls -R OUT

OUT:
ecDNA_prediction_DATA_0.95  ecDNA_summary_DATA_0.95

OUT/ecDNA_prediction_DATA_0.95:
LC675_AAACGAAAGGGTTCTT.txt  LC675_GTTACGAGTATGGGTG.txt
...

OUT/ecDNA_summary_DATA_0.95:
DATA_cnv.txt         DATA_pred.txt
DATA_count_freq.txt  DATA_ratio.txt
DATA_gini.txt
```
Seek `DATA_pred.txt` for the final prediction result.

Seek `DATA_count_freq.txt` for the contact frequency information summary.

### 2. Manhattan plot

Note: please run `CMPlot.R` in the interactive mode for flexibility. modify variable assignment `count_freq_file` to read count frequency file.

## Appendix

### 1. Directory hierarchy requirements

**Minimal requirements:**
```
> ls -R DATA

DATA:
LC675_AAACGAAAGGGTTCTT  LC675_CAGCCTTAGCTGATTC
...

DATA/LC675_AAACGAAAGGGTTCTT:
1000000.CNV.bedGraph                       matrix.mtx
```

### 2. File format requirements

#### 2.1 Single-cell copy number variation data
```
> head 1000000.CNV.bedGraph

chr1	0	1000000	0.0
chr1	1000000	2000000	1.4897193720663133
chr1	2000000	3000000	1.267909097945521
chr1	3000000	4000000	1.3012399754809931
chr1	4000000	5000000	1.2818908405862612
chr1	5000000	6000000	0.5364818463718226
chr1	6000000	7000000	0.8262895800095347
chr1	7000000	8000000	1.3324624442235051
chr1	8000000	9000000	0.9503856707800673
chr1	9000000	10000000	1.4120115495979793
```

#### 2.1 Single-cell copy Hi-C contact matrices
```
> head matrix.mtx

chrom1	start1	end1	chrom2	start2	end2	count
chr1	1000000	2000000	chr1	1000000	2000000	41
chr1	1000000	2000000	chr1	9000000	10000000	1
chr1	1000000	2000000	chr1	18000000	19000000	1
chr1	1000000	2000000	chr1	19000000	20000000	1
chr1	1000000	2000000	chr6	148000000	149000000	1
chr1	1000000	2000000	chr9	20000000	21000000	1
chr1	1000000	2000000	chr19	24000000	25000000	1
chr1	1000000	2000000	chr20	3000000	4000000	1
chr1	1000000	2000000	chr20	61000000	62000000	1
```

## References

* O. Tange (2011): GNU Parallel - The Command-Line Power Tool, The USENIX Magazine, February 2011:42-47.
* Yin, L. et al. rMVP: A Memory-efficient, Visualization-enhanced, and Parallel-accelerated tool for Genome-Wide Association Study, Genomics, Proteomics & Bioinformatics (2021), doi: 10.1016/j.gpb.2020.10.007.
* The pandas development team (2020). pandas-dev/pandas: Pandas (Version latest). doi:10.5281/zenodo.3509134.
* Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2.
* mckib2. (n.d.). Mckib2/Pygini: Compute the gini index. GitHub. https://github.com/mckib2/pygini

### Contact Us

For any question, contact Ming Hu (hum@ccf.org), Jiachen Sun (jxs2269@case.edu), or submit an issue on GitHub.
