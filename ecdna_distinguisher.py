import os
import random
import sys
from functools import partial

import numpy as np
import pandas as pd
import pygini
import torch
import torch.multiprocessing
import torch.nn as nn
import torch.nn.functional as F
from scipy.sparse import coo_matrix
from torch.multiprocessing import Pool

os.chdir(os.path.dirname(os.path.realpath(__file__)))

input_dir = sys.argv[1].rstrip('/')
output_dir = sys.argv[2].rstrip('/')
num_processes = int(sys.argv[3])
sample_name = input_dir.split("/")[-1]

chr_index = pd.read_csv("index.csv", header=0)

wsize = 5
msize = 3044

chr_size = [249, 243, 199, 191, 182, 171, 160, 146, 139, 134, 136, 134, 115, 108, 102, 91, 84, 81, 59, 65, 47, 51, 157]

random.seed(42)
torch.manual_seed(42)


def read_mtx(path, index):
    contact_mat = pd.read_csv(path, sep="\t", header=0)
    contact_mat = contact_mat[(contact_mat["chrom1"] != "chrY") & (contact_mat["chrom2"] != "chrY")]
    contact_mat = pd.concat([contact_mat["chrom1"] + "_" + contact_mat["start1"].astype(str),
                             contact_mat["chrom2"] + "_" + contact_mat["start2"].astype(str), contact_mat["count"]],
                            axis=1)
    contact_mat.columns = ["loc1", "loc2", "count"]
    contact_mat = pd.merge(contact_mat, index, left_on="loc1", right_on="loc", how="left").drop(["loc", "loc1"], axis=1)
    contact_mat = contact_mat.rename(columns={contact_mat.columns[2]: "index1"})
    contact_mat = pd.merge(contact_mat, index, left_on="loc2", right_on="loc", how="left").drop(["loc", "loc2"], axis=1)
    contact_mat = contact_mat.rename(columns={contact_mat.columns[2]: "index2"})
    contact_mat = coo_matrix((contact_mat['count'], (contact_mat['index1'], contact_mat['index2'])),
                             shape=(index.shape[0], index.shape[0]))
    contact_mat = contact_mat.toarray()

    return contact_mat


def slide(matrix, window_0, window_1):
    half_window_0 = int((window_0 - 1) / 2)
    half_window_1 = int((window_1 - 1) / 2)

    centers = []
    tensors = []

    for i in range(matrix.shape[0] - window_1 + 1):
        center = i + half_window_1
        lower = center - half_window_0
        upper = center + half_window_0 + 1

        tensors.append(matrix[lower:upper, :])
        centers.append(center)

    test_tensor = torch.stack(tensors, dim=0)
    centers = np.array(centers)

    return test_tensor, centers


def transform_file(matrix, chr_index, sample_name, output_dir, type="ecDNA"):
    matrix = matrix.astype(int)

    if type == "ecDNA":
        matrix = matrix.replace(2, 0)
        print("Summarizing ecDNA...")
    elif type == "HSR":
        matrix = matrix.replace(1, 0)
        matrix = matrix.replace(2, 1)
        print("Summarizing HSR...")
    elif type == "all":
        matrix = matrix.replace(2, 1)
        print("Summarizing all...")
    else:
        print("Invalid type.")

    freq = pd.DataFrame({"count": matrix.sum(axis=1), "freq": matrix.sum(axis=1) / matrix.shape[1]})

    freq = pd.concat([chr_index, freq], axis=1)

    freq[["chr", "start"]] = freq["loc"].str.split("_", expand=True)
    freq["end"] = freq["start"].astype(int) + 1000000

    freq = freq.drop(["loc", "index"], axis=1)

    freq = freq[freq.columns[2:].tolist() + freq.columns[:2].tolist()]

    freq.to_csv(output_dir + "/" + sample_name + "_summary_" + type + ".txt", index=False, sep="\t")


def summarize(results, chr_index, sample_name, output_dir, wsize):
    record = pd.concat([result for result in results if result is not None], axis=1)

    empty_rows = pd.DataFrame(0, columns=record.columns, index=range(int((wsize - 1) / 2)))
    record = pd.concat([empty_rows, record, empty_rows]).reset_index(drop=True)

    record.to_csv(output_dir + "/" + sample_name + "_bin_barcode_matrix" + ".txt", index=False, sep="\t")

    transform_file(record, chr_index, sample_name, output_dir, "ecDNA")
    transform_file(record, chr_index, sample_name, output_dir, "HSR")
    transform_file(record, chr_index, sample_name, output_dir, "all")


def process_file(model, chr_index, wsize, path):
    name = path.split("/")[-1]

    print("Processing",name)

    try:
        mat = read_mtx(path + "/matrix.mtx", chr_index)
    except:
        print("Error reading mtx file.")
        return

    if not mat.shape == (msize, msize):
        print("Matrix shape mismatch.")
        return

    test_tensor, centers = slide(torch.from_numpy(mat), wsize, wsize)

    with torch.no_grad():
        test_tensor = test_tensor.unsqueeze(1).float()

        test_pred = torch.argmax(model(test_tensor, centers), dim=1).detach().numpy()

        result = pd.concat([pd.DataFrame({name: test_pred}, index=range(2, len(test_pred) + 2))], axis=1)

        return result


class CNN(nn.Module):
    def __init__(self):
        super(CNN, self).__init__()

        self.conv_f2_l1 = nn.Conv2d(1, 8, kernel_size=(5, 45), stride=(1, 11), padding=(0, 2))
        self.conv_f2_l2 = nn.Conv2d(8, 16, kernel_size=(1, 45), stride=(1, 4), padding=(0, 0))

        self.bn_f2_l1 = nn.BatchNorm2d(8)
        self.bn_f2_l2 = nn.BatchNorm2d(16)
        self.bn_fc_l1 = nn.BatchNorm1d(64)

        self.pool_f2_l1 = nn.MaxPool2d(kernel_size=(1, 2), stride=2, padding=0)
        self.pool_f2_l2 = nn.MaxPool2d(kernel_size=(1, 2), stride=2, padding=0)

        self.fc_l1 = nn.Linear(5 * 5 + 12 * 16 + 5 + 1, 64)
        self.fc_l2 = nn.Linear(64, 3)

        self.dropout_f2 = nn.Dropout(0.5)
        self.dropout_fc = nn.Dropout(0.5)

        self.softmax = nn.Softmax(dim=1)

    def forward(self, x2, coord):

        start_indices = coord - int((5 + 1) / 2) + 1
        end_indices = coord + int((5 + 1) / 2)

        x2_chr = torch.split(x2[..., 2, :], chr_size, dim=-1)
        x2_c = torch.stack([torch.squeeze(torch.sum(st, dim=-1), dim=-1) for st in x2_chr], dim=-1)

        gini = []
        for i in range(x2_c.shape[0]):
            gini.append(torch.tensor(pygini.gini(x2_c[i, :].numpy()), dtype=torch.float32))
        gini = torch.unsqueeze(torch.stack(gini, dim=0), dim=-1)

        d4_slices = []

        for i in range(x2.size(0)):
            d4_slices.append(x2[i, ..., start_indices[i]:end_indices[i]])

        x1 = (torch.stack(d4_slices, dim=0) > 0).float()

        x2 = (x2 > 0).float()

        x2_rm = F.normalize(torch.squeeze(torch.sum(x2, dim=-1), dim=1), dim=-1, p=1)

        x2 = self.pool_f2_l1(F.relu(self.bn_f2_l1(self.conv_f2_l1(x2))))
        x2 = self.dropout_f2(x2)
        x2 = self.pool_f2_l2(F.relu(self.bn_f2_l2(self.conv_f2_l2(x2))))

        x1 = x1.view(-1, 5 * 5)
        x2 = x2.view(-1, 12 * 16)

        x = torch.cat((x1, x2, x2_rm, gini), dim=1)

        x = self.dropout_fc(F.gelu(self.bn_fc_l1(self.fc_l1(x))))

        x = self.softmax(self.fc_l2(x))

        return x


model = CNN()
model.load_state_dict(torch.load("model_dec20_2_multi.pt"))
model.eval()

dir0 = sorted([input_dir + "/" + d for d in os.listdir(input_dir) if ".DS_Store" not in d])

if __name__ == '__main__':
    print(f'Found {len(dir0)} cells...')

    with Pool(processes=num_processes) as pool:
        func = partial(process_file, model, chr_index, wsize)
        results = pool.map(func, dir0)

    summarize(results, chr_index, sample_name, output_dir, wsize)

    print("Sample", sample_name, "summarized.")
