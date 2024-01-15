# ecDNAcaller

**Last updated: Jan. 12, 2024**

## Prerequisites

### Packages

* **GNU Parallel**
* **Python**
    * **pygini**
    * **pandas**
    * **numpy**
    * **PyTorch** (for deep learning-based model only)
    * **scipy** (for deep learning-based model only)
* **R**
    * **CMplot**
* (other packages required as dependencies)

### Data

* **Copy number variation data** (for logistic regression-based model only)
* **Single-cell Hi-C contact matrices**

Note: see **File format requirements** for examples. File names are currently hard-coded in `ecDNAcaller` script. Modify
the section below if needed:

```
####################
cnv_name="1000000.CNV.bedGraph"
mat_name="matrix.mtx"
lm_dir=$script_dir"/coef_model_brain.txt"
####################
```

## Installation

We recommend to run `ecDNAcaller` and `ecDNAcaller_deep.py` in a conda/mamba environment.

```bash
conda create -n ecDNAcaller python parallel pandas numpy
pip install pygini
```

Then clone this repository to a local directory.

## Test Run

1. Download the example data, LC675 brain tumor dataset, and extract to a local directory.

2. Run command below (modify `-t` of `ecDNAcaller` and the last parameter of `ecDNAcaller_deep.py` according to the number of threads available):

```bash
sh ecDNAcaller -i example_data -o example_out_1 -p 0.95 -t 32

python ecDNAcaller_deep.py example_data example_out_2 32
```

3. Find the `example_data_count_freq.txt` in directory `example_out_1/ecDNA_summary_example_data_0.95` 

4. `example_data_summary_all.txt`, `example_data_summary_ecDNA.txt` and `example_data_summary_HSR.txt` in directory `example_out_2`.

5. Modify the following line with the actual path of files above in `CMPlot.R`:

```R
count_freq_file = "example_out/ecDNA_summary_example_data_0.95/example_data_count_freq.txt"
```

6. Then run `CMPlot.R` in interactive mode and generate plots.

## Usage

### 1. ecDNAcaller

```bash
sh ecDNAcaller -i <INPUT_PATH> -o <OUTPUT_PATH> -p <PROBABILITY_THRESHOLD> -t <THREADS>

sh ecDNAcaller -i <INPUT_PATH> -o <OUTPUT_PATH> -p <PROBABILITY_THRESHOLD> -t <THREADS> -s <BOOLEAN_OPTION>
```

**Arguments:**

`-i` input directory. See **Directory hierarchy requirements** for details.

`-o` output directory. See **Directory hierarchy requirements** for details.

`-p` probability cutoff threshold. 0.95 is recommended by default.

`-t` number of threads. Must be an integer.

`-s` (optional) summary-only mode. Use `-s true` to skip redundant cell processing (for example, if you need to change
the probability cutoff threshold and re-summarize).

Note: current version will ignore all interactions that involve chromosome Y (`'chrY'`).

**Output example:**

```
> ls -R example_out

example_out:
ecDNA_prediction_example_data  ecDNA_summary_example_data_0.95

example_out/ecDNA_prediction_example_data:
LC675_AAACGAAAGGGTTCTT.txt  LC675_GTTACGAGTATGGGTG.txt
...

example_out/ecDNA_summary_example_data_0.95:
example_data_cnv.txt         example_data_pred.txt
example_data_count_freq.txt  example_data_ratio.txt
example_data_gini.txt
```

Seek `example_data_pred.txt` for the final prediction result.

Seek `example_data_count_freq.txt` for the contact frequency information summary.

### 2. ecDNAcaller_deep.py

We developed a deep learning-based method that utilizes combined ideas from image classification as well as logistic
regression model above to achieve a more specific and efficient ecDNA/HSR detection.

With the same directory hierarchy, this deep learning-based model only requires
`matrix.mtx` for each cell as input.

```bash
python ecDNAcaller_deep.py <INPUT_PATH> <OUTPUT_PATH> <NUM_PROCESSES>
```

The script will output three summary files to the designated directory for 1) ecDNA alone, 2) HSR alone and 3) all together (that can also be
processed by `CMPlot.R` to generate a Manhattan plot) and a
bin-barcode binary matrix file in which row indices represent 10Mb bins from chr1 to chrX and column names are cell barcodes.

In the matrix, 0 represents *None*, 1 represents *ecDNA* and 2 represents *HSR*. Due to the algorithm design, the first 2 bins
of chromosome 1 and the last 2 bins of chromosome X will be padded with 0s.

The model currently runs on CPU only to allow for multiprocessing. Processing speed is at about 1.2 seconds per cell on
a single CPU core (Apple M1 Pro), which is approximately 5 times faster than the logistic regression model.

### 3. Manhattan plot

Note: please run `CMPlot.R` in interactive mode for flexibility. modify variable assignment `count_freq_file` to read
count frequency file.

## Appendix

### 1. Directory hierarchy requirements

**Minimal requirements:**

```
> ls -R example_data

example_data:
LC675_AAACGAAAGGGTTCTT  LC675_CAGCCTTAGCTGATTC
...

example_data/LC675_AAACGAAAGGGTTCTT:
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

Note: `_process.py` will add column names as:

```
'chr', 'start', 'end', 'cnv'
```

Thus, the columns must correspond.

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
* Yin, L. et al. rMVP: A Memory-efficient, Visualization-enhanced, and Parallel-accelerated tool for Genome-Wide
  Association Study, Genomics, Proteomics & Bioinformatics (2021), doi: 10.1016/j.gpb.2020.10.007.
* The pandas development team (2020). pandas-dev/pandas: Pandas (Version latest). doi:10.5281/zenodo.3509134.
* Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI:
  10.1038/s41586-020-2649-2.
* mckib2. (n.d.). Mckib2/Pygini: Compute the gini index. GitHub. https://github.com/mckib2/pygini

#### For deep learning-based model:

* Adam Paszke, Sam Gross, Francisco Massa, Adam Lerer, James Bradbury, Gregory Chanan, Trevor Killeen, Zeming Lin,
  Natalia Gimelshein, Luca Antiga, Alban Desmaison, Andreas Köpf, Edward Yang, Zach DeVito, Martin Raison, Alykhan
  Tejani, Sasank Chilamkurthy, Benoit Steiner, Lu Fang, Junjie Bai, and Soumith Chintala. 2019. PyTorch: an imperative
  style, high-performance deep learning library. Proceedings of the 33rd International Conference on Neural Information
  Processing Systems. Curran Associates Inc., Red Hook, NY, USA, Article 721, 8026–8037.
* Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski,
  Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod
  Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng,
  Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R
  Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020)
  SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.

## Contact Us

For any question, contact Ming Hu (hum@ccf.org), Jiachen Sun (jxs2269@case.edu), or submit an issue on GitHub.