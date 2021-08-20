import numpy as np
import pandas as pd
from scipy.stats import binom
import sys

# The code takes:
# Pileup_file
# Sample prefix

prefix = sys.argv[2] #Sample Prefix
total_data_read_file = sys.argv[1] #Pileup_file

total_data = pd.read_csv("%s" % total_data_read_file, header=None, sep="\t",
                         dtype={0: str, 1: int,
                                2: str, 3: int, 4: str, 5: str})

position = total_data[1]
ref_base = total_data[2]
coverage = total_data[3]
sequence = total_data[4]

nucs_matrix = np.zeros((6, len(total_data)), dtype=object)
G_histogram = np.empty((len(total_data), 5), dtype=object)
k = 0

for row_number in range(len(total_data)):
    flag = False
    c = len(sequence[row_number])
    base = ref_base[row_number]
    nucs_matrix[0][row_number] = position[row_number]
    nucs_matrix[1][row_number] = base
    for i in range(0, c):
        a = sequence[row_number][i]
        if flag == False:
            if a == "." or a == ",":
                if base == "A" or base == "a":
                    nucs_matrix[2][row_number] += 1
                elif base == "C" or base == "c":
                    nucs_matrix[3][row_number] += 1
                elif base == "G" or base == "g":
                    nucs_matrix[4][row_number] += 1
                elif base == "T" or base == "t":
                    nucs_matrix[5][row_number] += 1
            elif a == "A" or a == "a":
                nucs_matrix[2][row_number] += 1
            elif a == "C" or a == "c":
                nucs_matrix[3][row_number] += 1
            elif a == "G" or a == "g":
                nucs_matrix[4][row_number] += 1
            elif a == "T" or a == "t":
                nucs_matrix[5][row_number] += 1
            elif a == "+" or a == "-" or a == "^" or a == "$" or a == "*" or a == "<" or a == ">":
                flag = True
        elif a == "." or a == ",":
            flag = False
            if base == "A" or base == "a":
                nucs_matrix[2][row_number] += 1
            elif base == "C" or base == "c":
                nucs_matrix[3][row_number] += 1
            elif base == "G" or base == "g":
                nucs_matrix[4][row_number] += 1
            elif base == "T" or base == "t":
                nucs_matrix[5][row_number] += 1

    G_per = round(float(float(nucs_matrix[4][row_number]) / coverage[row_number]) * 100, 2)
    A_per = round(float(float(nucs_matrix[2][row_number]) / coverage[row_number]) * 100, 2)
    if base == "A" or base == "a" and G_per > 0:
            p = 0.001
            trials = round(float(float(G_per) * coverage[row_number]) / 100, 0)
            population = coverage[row_number]
            pvalue = round(1 - binom.cdf(trials, population, p, loc=0) + binom.pmf(trials, population, p, loc=0), 3)
            G_histogram[k][0] = position[row_number]
            G_histogram[k][1] = population
            G_histogram[k][2] = A_per
            G_histogram[k][3] = G_per
            G_histogram[k][4] = pvalue
            k += 1

G_Columns = ["Position", "Total Coverage", "% A Event", "% G Event", "pvalue"]
G_matrix = pd.DataFrame(G_histogram, columns=G_Columns).reset_index(
    drop=True).dropna()

nucs = ["Position", "Base", "A", "C", "G", "T"]
Motif_matrix = pd.DataFrame(nucs_matrix, index=nucs)

with pd.ExcelWriter("%s_Editing.xlsx" % prefix) as writer:
    Motif_matrix.to_excel(writer, sheet_name="Motif Matrix", header=False)
    G_matrix.to_excel(writer, sheet_name="A to G events", index=False)
