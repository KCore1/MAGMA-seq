import numpy as np
import random
import pandas as pd
import time
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
from scipy.optimize import minimize, fsolve
from scipy.special import gamma, erf
from scipy.stats import pearsonr, poisson
from regression import *

warnings.filterwarnings("ignore")


# Input parameters
sigma = 1.02
B = 350
kd_guesses = [50, 100, 200]
fmax_guesses = [40000, 30000, 25000]
conc_to_ignore = [0, 25]

bcs = [
    "TAGGACATCCGTAATTGTGG",
    "ACATGATTACTCATGGGTTC",
    "ATATGTAGCAATCTATCAAC",
    "CCGTCCTCATTTCCTTCTGT",
    "ATGCCTTCCTGGGTGGGTAC",
    "ACATATTCAAAGACGGCAGG",
    "TAGTAATTAAGTGATTCAAC",
    "ATGCACACATTTAAAGCTGT",
    "ACAGGCATGAATAATTATGG",
    "AAACCTACATTGCAAGGCTC",
    "ATTGCTGTGCGCGAATAATT",
    "TCGTCCATGATGGAAGCTTG",
    "ATGTATTTGTATCCTCGAGC",
    "TATGATGGCAACGTGTAATT",
    "CCATACATGCAGGAAGCAGC",
]

input_df = pd.read_csv(
    "/projects/brpe7306/fab-library-barcoding/Match_all/_combined.csv"
)
barcodes_to_test = pd.read_csv(
    "/projects/brpe7306/fab-library-barcoding/Match_all/_combined.csv"
)
barcodes_to_test = barcodes_to_test[barcodes_to_test["Barcode"].isin(bcs)]


# Run regression
output = run_reg_barcodes(
    input_df,
    barcodes_to_test,
    conc_to_ignore,
    sigma=sigma,
    B=B,
    kd_guesses=kd_guesses,
    fmax_guesses=fmax_guesses,
    top25_only=True,
    highmid_only=False,
    plot=False,
    print_df=False,
)
print(output)

output.to_csv("test_bc_reg.csv", index=False)
