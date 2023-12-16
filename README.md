# ecDNAfinder

**Last updated: November 15, 2023**

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

* **Copy number variation data**
* **Single-cell Hi-C contact matrices**
* **Trained logistic regression model**

Note: see **File format requirements** for examples. File names are currently hard-coded in `ecDNAfinder` script. Modify the section below if needed:
```
####################
cnv_name="1000000.CNV.bedGraph"
mat_name="matrix.mtx"
lm_dir=$script_dir"/coef_model_brain.txt"
####################
```

## Installation

We recommend to run `ecDNAfinder` in a conda/mamba environment.
```bash
conda create -n ecDNAfinder python parallel pandas numpy
pip install pygini
```
Then clone this repository to a local directory.

## Test Run

1. Download the example data, LC675 brain tumor dataset, and extract to a local directory.

2. Run command below (modify `-t` according to the number of threads available):
```bash
sh ecDNAfinder -i example_data -o example_out -p 0.95 -t 32
```
3. Find the `example_data_count_freq.txt` in directory `example_out/ecDNA_summary_example_data_0.95`  and modify the following line with the actual path in `CMPlot.R`:
```R
count_freq_file = "example_out/ecDNA_summary_example_data_0.95/example_data_count_freq.txt"
```
4. Then run `CMPlot.R` in interactive mode and generate a plot. The Manhattan plot should look like this:
   ![](https://github.com/Enterprise-D/ecDNAfinder/blob/main/Rect_Manhtn.LC675.jpg)

## Usage

### 1. ecDNAfinder

```bash
sh ecDNAfinder -i <INPUT_PATH> -o <OUTPUT_PATH> -p <PROBABILITY_THRESHOLD> -t <THREADS>

sh ecDNAfinder -i <INPUT_PATH> -o <OUTPUT_PATH> -p <PROBABILITY_THRESHOLD> -t <THREADS> -s <BOOLEAN_OPTION>
```

**Arguments:**

`-i` input directory. See **Directory hierarchy requirements** for details.

`-o` output directory. See **Directory hierarchy requirements** for details.

`-p` probability cutoff threshold. 0.95 is recommended by default.

`-t` number of threads. Must be an integer.

`-s` (optional) summary-only mode. Use `-s true` to skip redundant cell processing (for example, if you need to change the probability cutoff threshold and re-summarize).

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

### 2. Manhattan plot

Note: please run `CMPlot.R` in interactive mode for flexibility. modify variable assignment `count_freq_file` to read count frequency file.

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

## Alternative Method

We recently developed a neural network-based method to achieve a more specific and efficient ecDNA detection.

Before running the script, please have **PyTorch** and **scipy** installed in addition to the packages listed above.

```bash
python ecdna_detect_nn.py <INPUT_PATH> <OUTPUT_PATH> <PROB_CUTOFF> <NUM_PROCESSES>
```
The script will output a file which can also be processed by `CMPlot.R` to generate a Manhattan plot.

The model currently runs on CPU that allows for multiprocessing.

### Performance Comparison

The performance of prediction measured by ROC curve is significantly improved by the neural network model
on both our validation dataset (10% of our LC499/LC500 brain tumor dataset) and a separate dataset that represents colon cancer with distinct Hi-C pattern
that is hard to distinguish from ecDNA:

* **for validation dataset**:

![ROC1_NN.png](images%2FROC1_NN.png)

* **for LC676/LC677 dataset**:

![ROC2_NN.png](images%2FROC2_NN.png)

### Probability Cutoff Selection

Probability cutoff selection is a trade-off between sensitivity and specificity. 
For logistic regression, we chose 0.95 for minimized false positive rate. 
However, our current neural network model tends to reject more uncertain ecDNA candidates 
and pushes their predicted probability to a small number. Below is a comparison of the predicted probability on
a single cell (LC500_ACTAGGTGTTACCCAA) at different chromosomal bins (upper: logistic regression; lower: neural network):

* **linear model**:

![LC500_ACTAGGTGTTACCCAA_LM.png](images%2FLC500_ACTAGGTGTTACCCAA_LM.png)

* **neural model**:

![LC500_ACTAGGTGTTACCCAA_NN.png](images%2FLC500_ACTAGGTGTTACCCAA_NN.png)

On our validation dataset, we tested how different probability cutoffs affect prediction of positive and negative samples.
Figure below shows the difference between the true positive rate on positive samples and the false positive rate on negative samples.

![Delta.png](images%2FDelta.png)

Ideally this value should be maximized. However, choosing 0.05 as the cutoff will result in more false positives,
even though the true positive rate is at about 0.98 and the false positive rate is at about 0.02. 
Thus, we generally recommend a cutoff of 0.1 to 0.3 for the neural network model, depending on the desired sensitivity.
Figures below demonstrates effect on the Manhattan plot by choosing different cutoff using a mixed positive-negative dataset:

![Rect_Manhtn.0.05.jpg](images%2FRect_Manhtn.0.05.jpg)

![Rect_Manhtn.0.10.jpg](images%2FRect_Manhtn.0.10.jpg)

![Rect_Manhtn.0.15.jpg](images%2FRect_Manhtn.0.15.jpg)

![Rect_Manhtn.0.20.jpg](images%2FRect_Manhtn.0.20.jpg)

![Rect_Manhtn.0.30.jpg](images%2FRect_Manhtn.0.30.jpg)

![Rect_Manhtn.0.40.jpg](images%2FRect_Manhtn.0.40.jpg)

![Rect_Manhtn.0.50.jpg](images%2FRect_Manhtn.0.50.jpg)

## References

* O. Tange (2011): GNU Parallel - The Command-Line Power Tool, The USENIX Magazine, February 2011:42-47.
* Yin, L. et al. rMVP: A Memory-efficient, Visualization-enhanced, and Parallel-accelerated tool for Genome-Wide Association Study, Genomics, Proteomics & Bioinformatics (2021), doi: 10.1016/j.gpb.2020.10.007.
* The pandas development team (2020). pandas-dev/pandas: Pandas (Version latest). doi:10.5281/zenodo.3509134.
* Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2.
* mckib2. (n.d.). Mckib2/Pygini: Compute the gini index. GitHub. https://github.com/mckib2/pygini

#### For alternative method:

* Adam Paszke, Sam Gross, Francisco Massa, Adam Lerer, James Bradbury, Gregory Chanan, Trevor Killeen, Zeming Lin, Natalia Gimelshein, Luca Antiga, Alban Desmaison, Andreas Köpf, Edward Yang, Zach DeVito, Martin Raison, Alykhan Tejani, Sasank Chilamkurthy, Benoit Steiner, Lu Fang, Junjie Bai, and Soumith Chintala. 2019. PyTorch: an imperative style, high-performance deep learning library. Proceedings of the 33rd International Conference on Neural Information Processing Systems. Curran Associates Inc., Red Hook, NY, USA, Article 721, 8026–8037.
* Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.

## Contact Us

For any question, contact Ming Hu (hum@ccf.org), Jiachen Sun (jxs2269@case.edu), or submit an issue on GitHub.