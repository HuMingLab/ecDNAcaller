import os
import sys

import pandas as pd

input_dir = sys.argv[1]
output_dir = sys.argv[2]
prob_cutoff = sys.argv[3]

if input_dir.endswith('/'):
    input_dir = input_dir.rstrip('/')

if output_dir.endswith('/'):
    output_dir = output_dir.rstrip('/')

sample_name = input_dir.split('/')[-1]
prediction_dir = output_dir + "/" + "ecDNA_prediction_" + sample_name + '_' + prob_cutoff
summary_dir = output_dir + "/" + "ecDNA_summary_" + sample_name + '_' + prob_cutoff

# %%

prediction_file_list = os.listdir(prediction_dir)

cell_file_list = os.listdir(input_dir)
cell_file_list = [file for file in cell_file_list if not file.lower().endswith('.ds_store')]  # for macOS compatibility

# %%

a = pd.read_table(prediction_dir + "/" + prediction_file_list[0], header=0)
coord = a.iloc[:, 0:3]

# %%

m_cnv = pd.DataFrame(0, index=coord.index, columns=prediction_file_list)
m_ratio = m_cnv.copy()
m_gini = m_cnv.copy()
m_pred = m_cnv.copy()

# %%

# Loop through ann and read data from files
for i, cell in enumerate(prediction_file_list):
    res = pd.read_table(prediction_dir + "/" + cell, header=0)

    m_cnv[cell] = res['cnv']
    m_ratio[cell] = res['inter.intra.log2ratio']
    m_gini[cell] = res['gini']
    m_pred[cell] = res['pred']

# %%


# Concatenate output dataframes with 'out'
final_cnv = pd.concat([coord, m_cnv], axis=1)
final_cnv.columns = final_cnv.columns.str.replace('.txt', '')

final_ratio = pd.concat([coord, m_ratio], axis=1)
final_ratio.columns = final_ratio.columns.str.replace('.txt', '')

final_gini = pd.concat([coord, m_gini], axis=1)
final_gini.columns = final_gini.columns.str.replace('.txt', '')

final_pred = pd.concat([coord, m_pred], axis=1)
final_pred.columns = final_pred.columns.str.replace('.txt', '')

# %%

if not os.path.exists(summary_dir):
    os.makedirs(summary_dir)

# Write results to files
final_cnv.to_csv(f'{summary_dir}/{sample_name}_cnv.txt', sep='\t', index=False)
final_ratio.to_csv(f'{summary_dir}/{sample_name}_ratio.txt', sep='\t', index=False)
final_gini.to_csv(f'{summary_dir}/{sample_name}_gini.txt', sep='\t', index=False)
final_pred.to_csv(f'{summary_dir}/{sample_name}_pred.txt', sep='\t', index=False)

# %%

binary_pred = final_pred.iloc[:, 3:].copy()
for i in range(len(binary_pred.columns)):
    binary_pred[binary_pred.columns[i]] = (binary_pred[binary_pred.columns[i]] > float(prob_cutoff)).astype(int)

count = binary_pred.sum(axis=1)
freq = count / len(binary_pred.columns)

final_count_freq = pd.concat([coord, count, freq], axis=1)
final_count_freq.columns = ['chr', 'start', 'end', 'count', 'freq']

final_count_freq.to_csv(f'{summary_dir}/{sample_name}_count_freq.txt', sep='\t', index=False)

print("Sample", sample_name, "summarized.")
