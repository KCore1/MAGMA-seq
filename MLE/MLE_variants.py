import numpy as np
import random
import pandas as pd
import time
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
from scipy.optimize import minimize, fsolve
from scipy.special import gamma, erf, erfinv
from scipy.stats import pearsonr, poisson
from MLE_functions import *

warnings.filterwarnings("ignore")


# Input parameters
sigma = 1.02
B = 350
kd_guesses = [50, 70, 100]
fmax_guesses = [40000, 30000, 20000]
conc_to_ignore = [0]
input_csv = "./MLE/input/mandi_titration_100nM-HAB_repA.csv"
output_csv = "./MLE/output/mandi_100nM-HAB_repA_MLE.csv"

# Make sure WT is first variant for ddg calculations
variants = ["WT", "R59K", "S92V", "A108T", "S109R", "G122S", "M158G", "L159V"]

input_df = pd.read_csv(input_csv)
# input_df = input_df[input_df["Variant"].isin(variants)]  # comment out to run all variants

# Run MLE
start = time.time()
output = run_mle(
    input_df,
    conc_to_ignore,
    sigma=sigma,
    B=B,
    kd_guesses=kd_guesses,
    fmax_guesses=fmax_guesses,
    variant_col="Variant",
    top25_only=True,
    highmid_only=False,
    plot=False,
    print_df=False,
    verbose=False,
    status=True,
)
end = time.time()
output.to_csv(output_csv, index=False)
print(output)
print("Finished in", round((end - start), 2), "seconds.")
