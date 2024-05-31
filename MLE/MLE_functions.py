import numpy as np
import random
import pandas as pd
import time
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
import math
from scipy.optimize import minimize, fsolve, leastsq, basinhopping
from scipy.special import gamma, erf, erfinv
from scipy.stats import pearsonr, poisson, f

warnings.filterwarnings("ignore")

custom_params = {
    "axes.spines.right": False,
    "axes.spines.top": False,
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "figure.figsize": (7, 5),
}
sns.set(rc=custom_params, style="ticks", font_scale=1.5)


def get_gene(input_str):
    gene = input_str.split(">")[0]
    return gene


def extract_mutations_from_str(input_str):
    mut_list = input_str.split("|")
    mut_list = [mut.split(":")[1] for mut in mut_list]
    mut_list = [mut.split(";") for mut in mut_list]
    mut_list = mut_list[0] + mut_list[1]
    mut_list = [mut.split("-")[0] for mut in mut_list]
    mut_list = ",".join([mut for mut in mut_list if mut != "WT"])
    if mut_list == "":
        return "WT"
    return mut_list


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


def Fbar_func(**k):
    """Fmax: variant dependent maximum fluorescence encoded in units of RFU."""
    Fbar = (k["bright"] * k["Fmax"] - k["B"]) * k["Ck"] / (k["Kd"] + k["Ck"]) + k["B"]
    return max(Fbar, k["B"])


def pijk_func(**k):
    """Calculates probability of cell of variant i being sorted into bin j at label concentration k."""
    pijk = (1 - k["error_rate"]) * (
        0.5
        * erf(
            (np.log(k["Hj"]) - np.log(k["Fbar"]) + 0.5 * (k["sigma"] ** 2))
            / (k["sigma"] * np.sqrt(2))
        )
        - 0.5
        * erf(
            (np.log(k["Gj"]) - np.log(k["Fbar"]) + 0.5 * (k["sigma"] ** 2))
            / (k["sigma"] * np.sqrt(2))
        )
    ) + k["error_rate"]
    return pijk


def pijk_exp(**k):
    """Calculates mu, the mean of poisson that represents expected read counts."""
    return k["sorting_frac"] * (k["rijk"] / k["rjk"]) / k["fi"]


def sigmapijk_func(**k):
    sigma_pijk = np.sqrt(
        k['extrinsic_error']**2 + k["pijk_exp"] ** 2 * (1 / k["rijk"] + 1 / k["rir"])
    )
    return sigma_pijk


def loglikelihood_func(**k):
    return ((k["pijk_exp"] - k["pijk"]) / k["sigmapijk"]) ** 2


def calc_Fbar(**k):
    """Calculate the average fluorescence value given variant Kd and maximum fluorescence (Fmax)."""
    Fbar_ik = [Fbar_func(Kd=k["Kd"], Fmax=k["Fmax"], B=k["B"], Ck=k["Ck"][i], bright=k["bright"][i]) for i in range(len(k["Ck"]))]
    return Fbar_ik


def calc_pijk(**k):
    """Calculates the probability that a given variant will be sorted into a bin with fluorescence bounds Hj and Gj."""
    assert len(k["Fbar"]) == len(k["Hj"]) == len(k["Gj"])
    pijk_array = [
        pijk_func(
            Fbar=k["Fbar"][i],
            Hj=k["Hj"][i],
            Gj=k["Gj"][i],
            sigma=k["sigma"],
            error_rate=k["error_rate"],
        )
        for i in range(len(k["Fbar"]))
    ]
    return pijk_array


def calc_pijk_exp(**k):
    pijk_exp_array = [
        pijk_exp(
            fi=k["fi"],
            rijk=k["rijk"][i],
            rjk=k["rjk"][i],
            sorting_frac=k["sorting_frac"][i],
        )
        for i in range(len(k["rijk"]))
    ]
    return pijk_exp_array


def calc_sigmapijk(**k):
    sigma_pijk_array = [
        sigmapijk_func(pijk_exp=k["pijk_exp"][i], rijk=k["rijk"][i], rir=k["rir"], extrinsic_error=k['extrinsic_error'])
        for i in range(len(k["pijk_exp"]))
    ]
    return sigma_pijk_array


def calc_likelihood(**k):
    """Calculates the likelihood term for all variants and stores in array."""
    l_array = [
        loglikelihood_func(
            pijk=k["pijk"][i], pijk_exp=k["pijk_exp"][i], sigmapijk=k["sigmapijk"][i]
        )
        for i in range(len(k["pijk"]))
    ]
    return l_array


def sum_likelihood(solver_vars, B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, ext_error, bright):
    """Returns the negative of the summation of the log likelihood terms for each variant.
    Solver minimizes the result by changing solver_vars, an array of the sequential variant Kds,
    maximum fluorescence values, and underlying reference counts."""
    Kd_solver = solver_vars[0]
    Fmax_solver = solver_vars[1]
    Fbar_ik = calc_Fbar(Kd=Kd_solver, Fmax=Fmax_solver, B=B, Ck=Ck, bright=bright)
    pijk = calc_pijk(Fbar=Fbar_ik, Hj=Hj, Gj=Gj, sigma=sigma, error_rate=0)
    pijk_exp = calc_pijk_exp(fi=fi, rijk=rijk, rjk=rjk, sorting_frac=sorting_frac)
    sigma_pijk = calc_sigmapijk(pijk_exp=pijk_exp, rijk=rijk, rir=rir, extrinsic_error=ext_error)
    l_array = calc_likelihood(pijk=pijk, pijk_exp=pijk_exp, sigmapijk=sigma_pijk)
    return sum(l_array)


def sum_likelihood_fmax_only(Fmax, Kd, B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, ext_error, bright):
    Fbar_ik = calc_Fbar(Kd=Kd, Fmax=Fmax, B=B, Ck=Ck, bright=bright)
    pijk = calc_pijk(Fbar=Fbar_ik, Hj=Hj, Gj=Gj, sigma=sigma, error_rate=0)
    pijk_exp = calc_pijk_exp(fi=fi, rijk=rijk, rjk=rjk, sorting_frac=sorting_frac)
    sigma_pijk = calc_sigmapijk(pijk_exp=pijk_exp, rijk=rijk, rir=rir, extrinsic_error=ext_error)
    l_array = calc_likelihood(pijk=pijk, pijk_exp=pijk_exp, sigmapijk=sigma_pijk)
    return sum(l_array)


def plot_reconstruction(kd, Fmax, Ck, mf, variant):
    c = np.logspace(0, 4, 200)
    y = [hill(conc, Fmax, kd) for conc in c]  # kd*1000
    plt.plot(c, y)
    sns.scatterplot(x=Ck, y=mf, s=40)
    # textstr = '\n'.join((r'$F_{max}=%.0f \pm %.0f$' %(fit[0],perr[0]),
    #     r'$K_d=%.2f$' %(fit[1]))
    plt.xlabel("Concentration (nM)")
    plt.ylabel(r"$\bar{F_i}$")
    plt.xscale("log")
    plt.title(variant)
    plt.tight_layout()
    plt.savefig("plots/plot" + str(variant) + ".png")
    # print(variant, "fig saved")
    plt.clf()


def minimize_log_likelihood(initial_guess, bounds, B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, ext_error, bright):
    result = minimize(
        fun=sum_likelihood,
        x0=initial_guess,
        args=(B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, ext_error, bright),
        method="L-BFGS-B",
        bounds=bounds,
    )
    success = (
        True
        if (
            (bounds[0][0] + 2 < result["x"][0] < bounds[0][1] - 2)
            and not (initial_guess[0] - 1 < result["x"][0] < initial_guess[0] + 1)
            and (bounds[1][0] + 100 < result["x"][1] < bounds[1][1] - 100)
            and not (initial_guess[1] - 10 < result["x"][1] < initial_guess[1] + 10)
            and (result["success"])
        )
        else False
    )
    return result, success


def minimize_log_likelihood_fmax(initial_guess, Kd, B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, ext_error, bright):
    result = minimize(
        fun=sum_likelihood_fmax_only,
        x0=initial_guess,
        args=(Kd, B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, ext_error, bright),
        method="L-BFGS-B",
    )
    #     result = basinhopping(
    #         func=sum_likelihood_fmax_only,
    #         x0=initial_guess,
    #         minimizer_kwargs={
    #             "args": (
    #                 Kd,
    #                 B,
    #                 Ck,
    #                 Hj,
    #                 Gj,
    #                 sigma,
    #                 rijk,
    #                 rjk,
    #                 rir,
    #                 sorting_frac,
    #                 fi,
    #             ),
    #             "method": "L-BFGS-B",
    #         },
    #     )
    return result


def sweep_kd(val, B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, n, target, ext_error, bright):
    result = minimize_log_likelihood_fmax(20000, val, B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, ext_error, bright)
    return round(1 / (n - 1) * result["fun"] - target, 4)


def find_ci_roots(x0, B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, n, target, ext_error, bright):
    log_guesses = [
        (low, high)
        for low, high in zip([0.999, 0.8, 0.7, 0.6, 0.5], [1.01, 1.2, 1.3, 1.4, 1.5])
    ]
    for guess in log_guesses:
        ci_root_1 = fsolve(
            sweep_kd,
            x0 ** guess[0],
            args=(B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, n, target, ext_error, bright),
            xtol=0.001,
            full_output=True,
        )
        if ci_root_1[0] < x0 and ci_root_1[0] > 0 and ci_root_1[2] == 1:
            break
    for guess in log_guesses:
        ci_root_2 = fsolve(
            sweep_kd,
            x0 ** guess[1],
            args=(B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, n, target, ext_error, bright),
            xtol=0.001,
            full_output=True,
        )
        if ci_root_2[0] > x0 and ci_root_2[2] == 1:
            break
    ci_low, ci_high = round(ci_root_1[0][0], 2), round(ci_root_2[0][0], 2)
    if ci_root_1[2] != 1:
        ci_low = 0
    if ci_root_2[2] != 1:
        ci_high = 0
    if ci_low < 0:
        ci_low = 0
    if ci_high > 100_000:
        ci_high = 0
    return (ci_low, ci_high)


def find_ci_roots_step(x0, B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, n, target, ext_error, bright):
    if x0 < 0 or x0 > 9999:
        return (0, 0)
    guess_low = [round(x0 - i, 2) for i in range(1, 10000) if (x0 - i) > 0]
    for guess1 in guess_low:
        x1 = sweep_kd(
            guess1, B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, n, target, ext_error, bright
        )
        if x1 > 0:
            break
    guess_high = [round(x0 + i, 2) for i in range(1, 10000) if (x0 + i) < 10000]
    for guess2 in guess_high:
        x2 = sweep_kd(
            guess2, B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, n, target, ext_error, bright
        )
        if x2 > 0:
            break
    if guess1 < 1:
        guess1 = 0
    if guess2 > 9998:
        guess2 = 0
    return (guess1, guess2)


def objective_sweep_kd(x, *args):
    """Define the objective function to minimize for conversion
    from a root-finding problem to an optimization/minimization problem."""
    return np.square(sweep_kd(x, *args))


def find_ci_roots_optimizer(x0, B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, n, target, ext_error, bright):
    log_guesses = [
        (low, high)
        for low, high in zip([0.999, 0.8, 0.7, 0.6, 0.5], [1.01, 1.2, 1.3, 1.4, 1.5])
    ]

    minimizer_kwargs = {
        "method": "L-BFGS-B",
        "args": (B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, n, target, ext_error, bright),
    }

    # Perform the basinhopping optimization
    for guess in log_guesses:
        ci_root_1 = basinhopping(
            objective_sweep_kd,
            x0 ** guess[0],
            minimizer_kwargs=minimizer_kwargs,
            niter=100,
        )
        if (
            ci_root_1["x"][0] < x0
            and ci_root_1["x"][0] > 0
            and ci_root_1["success"] == True
        ):
            break
    for guess in log_guesses:
        ci_root_2 = basinhopping(
            objective_sweep_kd,
            x0 ** guess[1],
            minimizer_kwargs=minimizer_kwargs,
            niter=100,
        )
        if ci_root_2["x"][0] > x0 and ci_root_2["success"] == True:
            break
    ci_low, ci_high = round(ci_root_1["x"][0], 2), round(ci_root_2["x"][0], 2)
    if ci_root_1["success"] != 1:
        ci_low = 0
    if ci_root_2["success"] != 1:
        ci_high = 0
    if ci_low < 0:
        ci_low = 0
    if ci_high > 100_000:
        ci_high = 0
    return (ci_low, ci_high)


def sample_log_spaced(kd, num_points):
    log_spaced = np.logspace(1, 3.5, round(num_points), base=10, endpoint=False)
    return log_spaced


def plot_ci(kds_to_sample, x2, variant, target):
    plt.plot(kds_to_sample, x2)
    plt.xscale("log")
    plt.title(variant)
    plt.xlabel(r"$K_D$ (nM)")
    plt.ylabel(r"$\chi^2_R$")
    plt.plot(kds_to_sample, [target] * len(kds_to_sample))
    plt.tight_layout()

    if not os.path.exists("./plots"):
        os.makedirs("./plots")

    plt.savefig("./plots/CI" + str(variant) + ".png")
    # print(variant, "fig saved")
    plt.clf()


def run_mle(
    input_df,
    conc_to_ignore,
    sigma,
    B,
    ext_error,
    kd_guesses,
    fmax_guesses,
    variant_col,
    top25_only=False,
    highmid_only=False,
    mid_only=False,
    plot=False,
    print_df=False,
    verbose=False,
    status=True,
    brightness_adjust=False,
):
    (
        variants,
        bcs,
        kds,
        fmaxs,
        successes,
        logl,
        counts,
        concs,
        sf,
        pijks,
        modelijks,
        avg_counts,
        x2s,
        cil,
        cih,
        ns,
    ) = ([] for _ in range(16))
    total_vars = len(input_df[variant_col].unique())

    for i, bc in enumerate(input_df[variant_col].unique()):
        # Bin/concentration removal
        if top25_only:
            var_df = input_df[
                (input_df[variant_col] == bc)
                & (~input_df["Concentration"].isin(conc_to_ignore))
                & (input_df["Bin"] == "top25")
            ].copy()
        elif highmid_only:
            var_df = input_df[
                (input_df[variant_col] == bc)
                & (~input_df["Concentration"].isin(conc_to_ignore))
                & (
                    (input_df["Bin"] == "top25")
                    | ((input_df["Bin"] == "next25") & (input_df["Concentration"] > 99))
                )
            ].copy()
        elif mid_only:
            var_df = input_df[
                (input_df[variant_col] == bc)
                & (~input_df["Concentration"].isin(conc_to_ignore))
                & (input_df["Bin"] == "next25")
            ].copy()
        else:
            var_df = input_df[
                (input_df[variant_col] == bc)
                & (~input_df["Concentration"].isin(conc_to_ignore))
            ].copy()
            
        var_df = var_df.sort_values(['Concentration', 'Bin'])
        if status:
            print(
                str("".join(["█"] * int((i + 1) / total_vars * 100)))
                + " "
                + str(round((i + 1) / total_vars * 100, 1))
                + "% complete",
                end="\r",
            )

        if len(var_df.index) < 4:
            continue

        if print_df:
            print(var_df)

        # convert to list and extract to variables
        var_df = var_df.to_dict("list")
        variant = var_df["Variant"][0]
        if variant_col == "Barcode":
            rir = var_df["ref_barcode"][0]
        elif variant_col == "Variant":
            rir = var_df["ref_variant"][0] # NOTE: This is changed from "Fi"

        rijk = var_df["Rijk"]  # reads of variant i in bin j at labeling concentration k
        rjk = var_df["Rjk"]  # reads of bin j at labeling concentration k
        Hj = var_df["Hj"]  # fluorescence upper bound
        Gj = var_df["Gj"]  # fluorescence lower bound
        Ck = var_df["Concentration"]  # labeling concentration
        nk = var_df["Nk"]
        njk = var_df["Njk"]
        sumrir = var_df["ref_total"][0]
        if brightness_adjust:
            bright = var_df['brightness_adjust']
        else:
            bright = [1]*len(var_df)
        fi = rir / sumrir
        sorting_frac = [njk / nk for njk, nk in zip(njk, nk)]
        n = len(rijk)
        ns.append(n)
        conf_limits = 1 + (1 / (n - 1) * f.ppf(0.95, 1, n))
        sf.append([round(i, 2) for i in sorting_frac])
        counts.append(rijk)
        avg_counts.append(round(np.mean(rijk), 2))
        concs.append(Ck)

        initial_guesses = [[k, f] for k, f in zip(kd_guesses, fmax_guesses)]
        bounds = [(0.01, 20000), (5000, 90000)]

        # Run minimization
        for guess in initial_guesses:
            result, success = minimize_log_likelihood(guess, bounds, B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, ext_error, bright)
            solution = [round(result["x"][0], 2), round(result["x"][1], 2), round(result["fun"], 2)]
            initial_guess = guess
            if success == True:
                break
#         kds_to_sample = sample_log_spaced(round(solution[0]), num_points=200)
#         ll = []
#         for val in kds_to_sample:
#             sum_l = sum_likelihood([val, solution[1]], B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, ext_error)
#             ll.append(round(-1*sum_l, 3))
#         logl_output = pd.DataFrame({"Kd": kds_to_sample, "LL": ll})
#         logl_output.to_csv(variant + "_loglikelihood.csv", index=False)
        if verbose:
            print(solution)
        chi_sq_red_min = round(1 / (n - 2) * solution[2], 4)
        x2s.append(chi_sq_red_min)
        target = chi_sq_red_min * conf_limits

        # plot fluorescence reconstruction and CI reduced chi squared
        if plot:
            kds_to_sample = sample_log_spaced(round(solution[0]), num_points=200)
            x2= []
            for val in kds_to_sample:
                result = minimize_log_likelihood_fmax(20000, val, B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, ext_error, bright)
                x2.append(round(1 / (n - 1) * result["fun"], 2))
            if variant_col == "Barcode":
                plot_ci(kds_to_sample, x2, bc + "_" + variant, target)
            if variant_col == "Variant":
                plot_ci(kds_to_sample, x2, variant, target)


        # Solve for 95% confidence intervals - second step function solver if fsolve fails
        ci_roots = find_ci_roots(
            solution[0],
            B,
            Ck,
            Hj,
            Gj,
            sigma,
            rijk,
            rjk,
            rir,
            sorting_frac,
            fi,
            n,
            target,
            ext_error,
            bright,
        )
        if (ci_roots[0] == 0 or ci_roots[1] == 0) and solution[0] < 9999:
            #             ci_roots = find_ci_roots_step(solution[0], B, Ck, Hj, Gj, sigma, rijk, rjk, rir, sorting_frac, fi, n, target, ext_error)
            ci_roots = find_ci_roots_optimizer(
                solution[0],
                B,
                Ck,
                Hj,
                Gj,
                sigma,
                rijk,
                rjk,
                rir,
                sorting_frac,
                fi,
                n,
                target,
                ext_error,
                bright,
            )
        cil.append(ci_roots[0])
        cih.append(ci_roots[1])
        if verbose:
            print("conf_limits: ", conf_limits)
            print(r"X^2_red_min: ", chi_sq_red_min)
            print("Target: ", target)
            print("Kd:", solution[0])
            print("95% CI: ", ci_roots)
            print("ddg CI: ", [kd_to_ddg(ci_r, 100) for ci_r in ci_roots])
        bcs.append(bc)
        variants.append(variant)
        kds.append(solution[0])
        fmaxs.append(solution[1])
        logl.append(solution[2])
        successes.append(success)
        er = [round(np.log2((ri / rj) / fi), 2) for ri, rj in zip(rijk, rjk)]
        
        mean_f = calc_mean_fluor(sigma, Gj, sorting_frac, er)
        if plot:
            plot_reconstruction(solution[0], solution[1], Ck, mean_f, variant)
        # test
        Fbar_ik = calc_Fbar(Kd=solution[0], Fmax=solution[1], B=B, Ck=Ck, bright=bright)
        pijk = calc_pijk(Fbar=Fbar_ik, Hj=Hj, Gj=Gj, sigma=sigma, error_rate=0)
        pijk_exp = calc_pijk_exp(fi=fi, rijk=rijk, rjk=rjk, sorting_frac=sorting_frac)
        modelijks.append([round(p, 4) for p in pijk])
        pijks.append([round(p, 4) for p in pijk_exp])
    # output data to csv
    output_df = pd.DataFrame(
        {
            variant_col: bcs,
            "Variant": variants,
            "Kd": kds,
            "Fmax": fmaxs,
            "Success": successes,
            "LL": logl,
            "Rijk": counts,
            "Concentrations": concs,
            "Pijk": pijks,
            "Modelijk": modelijks,
            "sorted fractions": sf,
            "Avg_counts": avg_counts,
            "X^2_red_min": x2s,
            "95% CI low": cil,
            "95% CI high": cih,
            "n": ns,
        }
    )

    # ensure that WT is processed first
    kd_wt = 1_000_000     #1mM reference state
    output_df["dg"] = kd_to_ddg(output_df["Kd"], kd_wt)
    output_df["Ab"] = output_df["Variant"].apply(get_gene)
    output_df = output_df.sort_values(by=["Variant"])
    return output_df
