"""Functions to combine scanned barcode files into a single file."""
import numpy as np
import pandas as pd
import re


def read_barcode_files(file_names, limit_csv_path, ref_freqs):
    """Read in barcode files and add all the relevant data for MLE so all the data can be added into one file.

    file_names: List of file names of the matched csvs. Each should be in the form conc_bin_matched.csv
    limit_csv_path: Parameter that specifies the csv where the bin limits can be found
    ref_freqs: Dictionary of reference population frequencies for each barcode

    Returns a list of dataframes to be concatenated."""
    dfs = []
    limit_df = pd.read_csv(limit_csv_path)
    for f in file_names:
        df = pd.read_csv(f)
        # Store the label concentrations, bins, etc. from the name of the file
        names = re.search(r".*/(?P<concentration>[0-9.]+nM)_(?P<bin>.+)_matched.csv", f)
        if names is None:
            df.loc[:, "Concentration"] = np.nan
            df.loc[:, "Bin"] = re.match(r".*/(?P<bin>.+)_matched.csv", f).group("bin")
            fi = [
                ref_freqs[mut] if mut in ref_freqs else np.nan
                for mut in df["Variant"].values
            ]
            f_barcode = [
                ref_freqs[bc] if bc in ref_freqs else np.nan
                for bc in df["Barcode"].values
            ]
            df.loc[:, "ref_variant"] = fi
            df.loc[:, "ref_barcode"] = f_barcode
            # remainder will be NaN for the limits, nk, njk NOTE: is this what we want?
        else:
            # get high and low limit from input csv for specific conc and bin
            # appending conc and bin etc. to the end of the df

            # catch if values are not specified in the limit csv
            try:
                hj = limit_df.loc[
                    (limit_df["Concentration"] == names["concentration"])
                    & (limit_df["Bin"] == names["bin"]),
                    "Hj",
                ].values[0]
                gj = limit_df.loc[
                    (limit_df["Concentration"] == names["concentration"])
                    & (limit_df["Bin"] == names["bin"]),
                    "Gj",
                ].values[0]
                nk = limit_df.loc[
                    (limit_df["Concentration"] == names["concentration"])
                    & (limit_df["Bin"] == names["bin"]),
                    "Nk",
                ].values[0]
                njk = limit_df.loc[
                    (limit_df["Concentration"] == names["concentration"])
                    & (limit_df["Bin"] == names["bin"]),
                    "Njk",
                ].values[0]
            except IndexError:
                print(
                    f"Error: Could not find limits for {names['concentration']} and {names['bin']} in {limit_csv_path}"
                )

            fi = [
                ref_freqs[mut] if mut in ref_freqs else np.nan
                for mut in df["Variant"].values
            ]
            f_barcode = [
                ref_freqs[bc] if bc in ref_freqs else np.nan
                for bc in df["Barcode"].values
            ]
            ref_total = ref_freqs["total"]
            df.loc[:, "Concentration"] = names[
                "concentration"
            ]  # Assuming micromolar units
            df.loc[:, "Bin"] = names["bin"]
            df.loc[:, "Hj"] = hj
            df.loc[:, "Gj"] = gj
            df.loc[:, "Nk"] = nk
            df.loc[:, "Njk"] = njk
            df.loc[:, "ref_variant"] = fi
            df.loc[:, "ref_barcode"] = f_barcode
            df.loc[:, "ref_total"] = ref_total
        dfs.append(df)
    return dfs
