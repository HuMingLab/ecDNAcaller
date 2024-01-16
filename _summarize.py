import os
import sys
import time

import pandas as pd

input_dir = sys.argv[1].rstrip('/')
output_dir = sys.argv[2].rstrip('/')
prob_cutoff = sys.argv[3]

sample_name = input_dir.split('/')[-1]
cache_dir = output_dir + "/" + sample_name + "_cache"
summary_dir = output_dir + "/" + sample_name + "_summary_" + prob_cutoff + "_prob"

print("Sample", sample_name, "| Summarizing ecDNA...")

# %%

cache_file_list = sorted(os.listdir(cache_dir))
cache_file_list = [file for file in cache_file_list if not file.lower().endswith('.ds_store')]

cell_file_list = sorted(os.listdir(input_dir))
cell_file_list = [file for file in cell_file_list if not file.lower().endswith('.ds_store')]  # for macOS compatibility

# %%

a = pd.read_table(cache_dir + "/" + cache_file_list[0], header=0)
coord = a.iloc[:, 0:3]

# %%

file_paths = [os.path.join(cache_dir, cell) for cell in cache_file_list]

dfs = [pd.read_table(file_path, header=0) for file_path in file_paths]

all_data = pd.concat(dfs, axis=1)

m_cnv = all_data[['cnv']]
m_ratio = all_data[['log2ratio']]
m_gini = all_data[['gini']]
m_pred = all_data[['pred']]

m_cnv.columns = cache_file_list
m_ratio.columns = cache_file_list
m_gini.columns = cache_file_list
m_pred.columns = cache_file_list

# %%

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

binary_pred = final_pred.iloc[:, 3:]
binary_pred = (binary_pred > float(prob_cutoff)).astype(int)

count = binary_pred.sum(axis=1)
freq = count / len(binary_pred.columns)

final_count_freq = pd.concat([coord, count, freq], axis=1)
final_count_freq.columns = ['chr', 'start', 'end', 'count', 'freq']

try:
    final_cnv.to_csv(f'{summary_dir}/{sample_name}_bin_barcode_matrix_cnv.txt', sep='\t', index=False)
    final_ratio.to_csv(f'{summary_dir}/{sample_name}_bin_barcode_matrix_ratio.txt', sep='\t', index=False)
    final_gini.to_csv(f'{summary_dir}/{sample_name}_bin_barcode_matrix_gini.txt', sep='\t', index=False)
    final_pred.to_csv(f'{summary_dir}/{sample_name}_bin_barcode_matrix_pred.txt', sep='\t', index=False)
    final_count_freq.to_csv(f'{summary_dir}/{sample_name}_count_freq_ecDNA.txt', sep='\t', index=False)

except KeyboardInterrupt:
    print("Sample", sample_name, "| File writing interrupted.")

    os.remove(f'{summary_dir}/{sample_name}_cnv.txt')
    os.remove(f'{summary_dir}/{sample_name}_ratio.txt')
    os.remove(f'{summary_dir}/{sample_name}_gini.txt')
    os.remove(f'{summary_dir}/{sample_name}_pred.txt')
    os.remove(f'{summary_dir}/{sample_name}_count_freq.txt')

    time.sleep(0.5)
    print("Sample", sample_name, "| Cache cleared.")

    exit(1)

print("Sample", sample_name, "| Summarized.")
