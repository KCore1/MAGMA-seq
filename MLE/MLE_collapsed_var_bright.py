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
ext_error = 0.02
kd_guesses = [50, 70, 100]
fmax_guesses = [40000, 30000, 20000]
conc_to_ignore = [0]

# input_csv = "/projects/brpe7306/fab-library-barcoding/Match_final_concat/_collapsed.csv"
# input_csv = "/projects/brpe7306/fab-library-barcoding/Match_final/_collapsed.csv"
# input_csv = "/projects/brpe7306/fab-library-barcoding/Match_COV2_final/_collapsed.csv"
input_csv = "/projects/brpe7306/fab-library-barcoding/MLE/input/rep1/HA_rep1_collapsed_bright_adjust.csv"

output_csv = "/projects/brpe7306/fab-library-barcoding/MLE/output/topbinHA_rep1_bright_test.csv"

input_df = pd.read_csv(input_csv)
variants = input_df['Variant'].unique()[:10]
input_df = input_df[input_df['Variant'].isin(variants)]

# Run MLE
start = time.time()
output = run_mle(
    input_df,
    conc_to_ignore,
    sigma=sigma,
    B=B,
    ext_error=ext_error,
    kd_guesses=kd_guesses,
    fmax_guesses=fmax_guesses,
    variant_col="Variant",
    top25_only=True,
    highmid_only=False,
    mid_only=False,
    plot=False,
    print_df=True,
    verbose=False,
    brightness_adjust=True,
)
end = time.time()
output.to_csv(output_csv, index=False)
print(output)
print("Finished in", round((end - start), 2), "seconds.")
