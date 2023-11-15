import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
# don't know why but solved the large-scale parallelization issue on HPC
# https://stackoverflow.com/questions/52026652/openblas-blas-thread-init-pthread-create-resource-temporarily-unavailable

import sys

import numpy as np
import pandas as pd
import pygini

# %%
cell_dir = sys.argv[1]
script_dir = sys.argv[2]
lm_dir = sys.argv[3]
cnv_name = sys.argv[4]
mat_name = sys.argv[5]
input_dir = sys.argv[6]
output_dir = sys.argv[7]

# %%

if cell_dir.endswith('/'):
    cell_dir = cell_dir.rstrip('/')

if script_dir.endswith('/'):
    script_dir = script_dir.rstrip('/')

if input_dir.endswith('/'):
    input_dir = input_dir.rstrip('/')

if output_dir.endswith('/'):
    output_dir = output_dir.rstrip('/')

# %%

cell_name = cell_dir.split('/')[-1]
sample_name = input_dir.split('/')[-1]

path_bed_graph = cell_dir + "/" + cnv_name
path_contact_matrix = cell_dir + "/" + mat_name

prediction_dir = output_dir + "/" + "ecDNA_prediction_" + sample_name
output_file_path = prediction_dir + '/' + cell_name + '.txt'

if (not os.path.exists(path_bed_graph)) or (not os.path.exists(path_contact_matrix)):
    print("Cell", cell_name, "error: one or more file(s) not exist. Skipped.")
    exit(1)

if os.path.exists(output_file_path):
    print("Cell", cell_name, "error: already processed. Skipped.")
    exit(3)

# %%
try:
    res = pd.read_table(path_bed_graph, header=None)
    mat = pd.read_table(path_contact_matrix)
except Exception as e:
    # Print or log the specific error message
    print("Cell", cell_name, "error: one or more file(s) not valid. Skipped.")
    # Raise a custom error
    exit(1)

if len(res) == 0 or len(mat) == 0:
    print("Cell", cell_name, "error: one or more file(s) not valid. Skipped.")
    exit(1)

# %%
num_chr = 23
name_chr = ['chr' + str(i) for i in range(1, 23)] + ['chrX']

res.columns = ['chr', 'start', 'end', 'cnv']

mat = mat[(mat['chrom1'] != 'chrY') & (mat['chrom2'] != 'chrY')]
res = res[res['chr'] != 'chrY']

res['num.intra.bin'] = 0
res['num.inter.bin'] = 0
# %%
rec = pd.DataFrame(np.zeros((len(res), num_chr), dtype=int), columns=name_chr)
# %%

# Any way to vectorize this?
for j in range(len(res)):
    overlap = mat[(mat['chrom1'] == res['chr'][j]) & (mat['start1'] == res['start'][j]) | (mat['chrom2'] == res['chr'][j]) & (mat['start2'] == res['start'][j])]
    intra = overlap[(overlap['chrom1'] == overlap['chrom2']) & (overlap['start1'] == res['start'][j])]
    inter = overlap[(overlap['chrom1'] != overlap['chrom2']) & (overlap['chrom1'] == res['chr'][j])]

    #inter = inter[(inter['chrom1'] != 'chrY') & (inter['chrom2'] != 'chrY')]

    res.at[j, 'num.intra.bin'] = len(intra)
    res.at[j, 'num.inter.bin'] = len(inter)

    for k in range(num_chr):
        rec.at[j, name_chr[k]] = np.sum(inter['chrom2'] == name_chr[k])

# %%
res['gini'] = np.nan

gini_values = np.apply_along_axis(lambda x: pygini.gini(x.astype("float64")), axis=1, arr=rec.values)
non_zero_sum_rows = rec.sum(axis=1) > 0
res.loc[non_zero_sum_rows, 'gini'] = gini_values[non_zero_sum_rows]
# %%
res['log2ratio'] = round(np.log2(res['num.inter.bin'] + 1) - np.log2(res['num.intra.bin'] + 1), 4)
# %%
linear_model = pd.read_csv(lm_dir, sep='\t', header=0)
# %%
input = res.iloc[:, [0, 1, 2, 3, 7, 6]]
input.columns = ['chr', 'start', 'end', 'cnv', 'ratio', 'gini']
input.loc[:, 'cnv'] = input['cnv'].apply(lambda value: min(value, 100))
# %%
res['eta'] = linear_model['Estimate'][0] + linear_model['Estimate'][1] * input['cnv'] + linear_model['Estimate'][2] * \
             input['ratio'] + linear_model['Estimate'][3] * input['gini']
res['pred'] = 1 / (1 + np.exp(-res['eta']))

res.loc[res['gini'].isna(), 'pred'] = 0

# %%
if not os.path.exists(prediction_dir):
    os.makedirs(prediction_dir)

res.to_csv(output_file_path, sep='\t', index=False, header=True, quoting=0)

print("Cell", cell_name, "processed. Exiting.")
