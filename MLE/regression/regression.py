import numpy as np
import random
import pandas as pd
import time
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
from scipy.optimize import curve_fit
from scipy.special import erf, erfinv
from scipy.stats import pearsonr, poisson

warnings.filterwarnings("ignore")

custom_params = {
    "axes.spines.right": False,
    "axes.spines.top": False,
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "figure.figsize": (8, 6),
}
sns.set(rc=custom_params, style="ticks", font_scale=2)


def kd_to_ddg(kd_mut, kd_ref):
    """Calculates the change in free energy from a reference to a mutant."""
    R = 0.001987  # gas constant (kcal/(mol*K))
    T = 310  # temperature (K)
    dg_ref = -R * T * np.log(kd_ref)
    dg_mut = -R * T * np.log(kd_mut)
    ddg = round(dg_ref - dg_mut, 3)
    return ddg


def MSE(x, y):
    """Calculates mean square error of arrays x and y"""
    se = [(x - y) ** 2 for x, y in zip(x, y)]
    mse = np.mean(se)
    return round(mse, 3)


def hill(x, fmax, kd):
    """Hill Equation"""
    n = 1
    return fmax * x**n / (kd**n + x**n)


def mean_fluor_func(sigma, Fg, sort_frac, ER):
    f_bar_i = np.exp(
        sigma**2 / 2
        + np.log(Fg)
        - sigma * np.sqrt(2) * erfinv(1 - sort_frac * 2 ** (ER + 1))
    )
    return round(f_bar_i, 2)


def calc_mean_fluor(sigma, Fg, sort_frac, ER):
    fbars = [mean_fluor_func(sigma, Fg[i], sort_frac[i], ER[i]) for i in range(len(Fg))]
    return fbars


def plot_reconstruction(kd, Fmax, Ck, mf, variant):
    Ck = [ck for ck in Ck]  # 1000*ck
    c = np.logspace(0, 5, 200)
    y = [hill(conc, Fmax, kd) for conc in c]  # kd*1000
    plt.plot(c, y)
    sns.scatterplot(x=Ck, y=mf)
    # textstr = '\n'.join((r'$F_{max}=%.0f \pm %.0f$' %(fit[0],perr[0]),
    #     r'$K_d=%.2f$' %(fit[1]))
    plt.xlabel("Concentration (pM)")
    plt.ylabel(r"$\bar{F_i}$")
    plt.xscale("log")
    plt.title(variant)
    plt.tight_layout()
    plt.savefig("plots/plot" + str(variant) + ".png")
    print(variant, "fig saved")
    plt.clf()


def run_reg_barcodes(
    input_df,
    barcodes_to_test,
    conc_to_ignore,
    sigma,
    B,
    kd_guesses,
    fmax_guesses,
    top25_only=False,
    highmid_only=False,
    plot=False,
    print_df=False,
):
    bcs, counts, concs, sf, ers, avg_counts, variants, mean_fs = ([] for _ in range(8))

    for i, bc in enumerate(barcodes_to_test["Barcode"].unique()):
        if top25_only:
            bc_df = input_df[
                (input_df["Barcode"] == bc)
                & (~input_df["Concentration"].isin(conc_to_ignore))
                & (input_df["Bin"] == "top25")
            ].copy()
        elif highmid_only:
            bc_df = input_df[
                (input_df["Barcode"] == bc)
                & (~input_df["Concentration"].isin(conc_to_ignore))
                & (
                    (input_df["Bin"] == "top25")
                    | ((input_df["Bin"] == "next25") & (input_df["Concentration"] > 99))
                )
            ].copy()
        else:
            bc_df = input_df[
                (input_df["Barcode"] == bc)
                & (~input_df["Concentration"].isin(conc_to_ignore))
            ].copy()

        if print_df:
            print(bc_df)

        try:
            # convert to list and extract to variables
            bc_df = bc_df.to_dict("list")
            rijk = bc_df[
                "Rijk"
            ]  # number of reads of variant i in bin j at labeling concentration k
            rjk = bc_df["Rjk"]  # number of reads of bin j at labeling concentration k
            Hj = bc_df["Hj"]  # fluorescence upper bound
            Gj = bc_df["Gj"]  # fluorescence lower bound
            Ck = bc_df["Concentration"]  # labeling concentration
            fi = bc_df["F_barcode"][0]  # reference frequency (float)
            variant = bc_df["Mutations"][0]
            sorting_frac = [njk / nk for njk, nk in zip(bc_df["Njk"], bc_df["Nk"])]

            initial_guesses = [[k, f] for k, f in zip(kd_guesses, fmax_guesses)][0]
            bounds = [(1, 20000), (5000, 90000)]

            sf.append([round(i, 2) for i in sorting_frac])
            counts.append(rijk)
            avg_counts.append(np.mean(rijk))
            concs.append(Ck)
            bcs.append(bc)
            variants.append(variant)
            er = [round(np.log2((ri / rj) / fi), 2) for ri, rj in zip(rijk, rjk)]
            ers.append(er)
            mean_f = calc_mean_fluor(sigma, Gj, sorting_frac, er)
            mean_fs.append(mean_f)

            # curve_fit
            popt, pcov = curve_fit(hill, Ck, mean_f, initial_guesses)
            print(variant)
            print(mean_f)
            print(popt)
            print(pcov)

            break
            plot_reconstruction(100, 10000, Ck, mean_f, variant)

        except:
            continue
    output_df = pd.DataFrame(
        {
            "Barcode": bcs,
            "Variant": variants,
            "Rijk": counts,
            "Concentrations": concs,
            "sorted fractions": sf,
            "ER": ers,
            "Avg_counts": avg_counts,
            "Mean_fluorescence": mean_fs,
        }
    )
    kd_wt = 100
    output_df["ddg"] = kd_to_ddg(100, kd_wt)
    return output_df
