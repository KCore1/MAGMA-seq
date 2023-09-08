import os
import numpy as np
import pandas as pd
import sys
import re
import time
from utilities.haplotype import open_by_extension
from utilities.merge import read_seqs, byte_array_to_str
from utilities.arguments import parse_args_and_read_config_match
import doctest
from utilities.read_barcode_files import read_barcode_files
import multiprocessing as mp

MIXED_BASES = {
    "A": ["A"],
    "T": ["T"],
    "C": ["C"],
    "G": ["G"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "M": ["A", "C"],
    "K": ["G", "T"],
    "S": ["C", "G"],
    "W": ["A", "T"],
    "H": ["A", "C", "T"],
    "B": ["C", "G", "T"],
    "V": ["A", "C", "G"],
    "D": ["A", "G", "T"],
    "N": ["A", "C", "G", "T"],
}


def barcode(s, template, start):
    """Extract barcode with regex, checking if the barcode matches the template

    s -- A read (string or byte array).
    template -- A barcode template (string).
    start -- Index of the first base of the barcode.

    Defining doctests:
    >>> s = "AGATCCAGC"
    >>> template = "ATC"
    >>> barcode_start = 2
    >>> barcode(s, template, barcode_start)
    'ATC'

    >>> s = "AGATACAGC"
    >>> barcode(s, template, barcode_start)

    >>> s = "AGATCCAGC"
    >>> template = "ATS"
    >>> barcode_start = 2
    >>> barcode(s, template, barcode_start)
    'ATC'

    """
    # NOTE: This is not the case any more with indexing -> There might be some case where the constant region occurs multiple times,
    # but I think that's really unlikely. If that happens, this will not identify the barcode properly.
    barcode = byte_array_to_str(s[start : start + len(template)])
    # barcode = re.search(r"{}(.{{{}}})".format(constant_reg, len(template)), s).group(1)
    for base, template_base in zip(barcode, template):
        if base not in MIXED_BASES[template_base]:
            return None
    return barcode


def unique_variants(df, d=dict()):
    """Get the unique variant frequencies.

    d: dictionary of variants and their counts (typically empty, but can be used to add to existing dictionary)
    df: dataframe of matched barcodes and variants/mutations
    TODO: Check if this is what we want, do we want this in the result? Or in the haplotyping code? Do we want below?
    """
    for var in df["Variant"].values:
        if var in d:
            d[var] += 1
        else:
            d[var] = 1
    total = sum(d.values())
    d = {k: v / total for k, v in d.items()}
    return d


def barcodes_per_variant(df, d=dict()):
    """Get the number of barcodes per variant/mutation.

    df: dataframe of matched barcodes and variants/mutations
    d: dictionary of variants and barcodes associated (typically empty, but can be used to add to existing dictionary)
    TODO: Check if this is what we want, do we want this in the result? Or in the haplotyping code? Do we want above?
    """
    for var in df["Variant"].values:
        if var in d:
            d[var].append(df[df["Variant"] == var]["Barcode"].values[0])
        else:
            d[var] = [df[df["Variant"] == var]["Barcode"].values[0]]
    d = {k: (v, len(v)) for k, v in d}
    return d


def match_barcodes(
    barcodes: pd.DataFrame, file, template, barcode_start, sample_name, debug=False
):
    """Match barcodes to reads and return a list of matches.

    barcodes -- A barcode variant map (DataFrame: Barcode, Variant, Freq, Count, P-value).
    file -- An open file handle of barcodes from sorted population (read from open file handle from gzip.open).
    template -- A barcode template (string).
    barcode_start -- The start of the barcode (int).

    Defining doctests:
    >>> barcodes = pd.DataFrame({'Barcode': ['ATG', 'ATC'], 'Variant': ['A', 'C'], 'Freq': [0.5, 0.5], 'Count': [1, 1], 'P-value': [0.5, 0.5]})
    >>> reads = ['GAATGCAB', 'GAATGCAC', 'GAATGCAT', 'GAATGCAG', 'GAATACGA', 'GAATACGC', 'GAATCCGT', 'GAATACGG']
    >>> template = 'ATS'
    >>> barcode_start = 2
    >>> match_barcodes(barcodes, reads, template, barcode_start)
    {'ATG': ['A', 'A', 'A', 'A'], 'ATC': ['C']}

    """
    num_unmatched = 0
    matched = dict()
    if type(file) == list:
        for seq in file:
            read_barcode = barcode(seq, template, barcode_start)
            if read_barcode in barcodes["Barcode"].values:
                if read_barcode in matched:
                    matched[read_barcode].append(
                        barcodes[barcodes["Barcode"] == read_barcode]["Variant"].values[
                            0
                        ]
                    )
                else:
                    matched[read_barcode] = [
                        barcodes[barcodes["Barcode"] == read_barcode]["Variant"].values[
                            0
                        ]
                    ]
    else:  # file is a file handle
        for id, seq, space, qual in read_seqs(file):
            read_barcode = barcode(seq, template, barcode_start)
            if read_barcode in barcodes["Barcode"].values:
                # print(read_barcode)
                if read_barcode in matched:
                    matched[read_barcode].append(
                        barcodes[barcodes["Barcode"] == read_barcode]["Variant"].values[
                            0
                        ]
                    )
                else:
                    matched[read_barcode] = [
                        barcodes[barcodes["Barcode"] == read_barcode]["Variant"].values[
                            0
                        ]
                    ]
            else:
                num_unmatched += 1
    if debug:
        total_matched_reads = sum([len(v) for v in matched.values()])
        print(
            f"Sample: {sample_name}\nNumber of matched barcodes = {len(matched)}\nPercent of matched reads = {round(100*total_matched_reads/(total_matched_reads+num_unmatched),2)}\n"
        )
    # this second parameter is rjk in MLE
    return matched, num_unmatched + total_matched_reads


def match_sample(sample, params, barcodes, ref_freqs, filenames):
    """Match barcodes for a sample."""
    print("Matching barcodes for sample {}...".format(sample[0]))
    with open_by_extension(
        os.path.join(params.fastq_file_dir, sample[1]), "rt"
    ) as file:
        matched, rjk = match_barcodes(
            barcodes, file, params.template, params.barcode_start, sample[0], True
        )  # debug is true currently
    lengths = [len(x) for x in matched.values()]
    for key in matched.keys():
        matched[key] = matched[key][0]
    matched = pd.DataFrame(
        {
            "Barcode": list(matched.keys()),
            "Variant": list(matched.values()),
            "Count": lengths,
            "Rjk": [rjk] * len(matched.keys()),
        }
    )
    # if samples are names properly, this will capture the concentration and the bin
    if "ref" in sample[0]:
        # ref_freqs/fi should be based on variant
        for _, mut in matched["Variant"].items():
            # ref_freqs[mut] = sum(matched[matched['Variant'] == mut]
            #                      ['Count'].values) / sum(matched['Count'].values)
            # change these to counts only
            ref_freqs[mut] = sum(matched[matched["Variant"] == mut]["Count"].values)
            ref_freqs["total"] = sum(matched["Count"].values)  # Get this into final csv
        # but we want barcode frequencies too for barcode enrichment ratio
        for _, bc in matched["Barcode"].items():
            # ref_freqs[bc] = sum(matched[matched['Barcode'] == bc]
            #                     ['Count'].values) / sum(matched['Count'].values)
            # change these to counts only
            ref_freqs[bc] = sum(matched[matched["Barcode"] == bc]["Count"].values)
        matched.to_csv(os.path.join(params.output_dir, "ref_matched.csv"), index=False)
    elif "conc" in sample[0]:
        conc = re.search(r"conc([0-9.]+nM)", sample[0]).group(1)
        bin = re.search(r".*bin(.+)", sample[0]).group(1)
        fn = os.path.join(params.output_dir, "{}_{}_matched.csv".format(conc, bin))
        matched.to_csv(
            os.path.join(params.output_dir, "{}_{}_matched.csv".format(conc, bin)),
            index=False,
        )
        filenames.append(fn)
    else:
        fn = os.path.join(params.output_dir, "{}_matched.csv".format(sample[0]))
        matched.to_csv(
            os.path.join(params.output_dir, "{}_matched.csv".format(sample[0])),
            index=False,
        )
        filenames.append(fn)


def matching(params, samples, barcodes):
    if params.use_multiprocessing:
        pool = mp.Pool(mp.cpu_count())
        filenames = mp.Manager().list()  # capturing filenames for combining
        ref_freqs = mp.Manager().dict()  # capturing reference population frequencies
        args = [(s, params, barcodes, ref_freqs, filenames) for s in samples.items()]
        pool.starmap_async(match_sample, args)
        pool.close()
        pool.join()

    else:
        filenames = []  # capturing filenames for combining
        ref_freqs = dict()  # capturing reference population frequencies

        for s in samples.items():
            match_sample(s, params, barcodes, ref_freqs, filenames)

    return filenames, ref_freqs


def output_collapsed_counts(df):
    """Output a collapsed dataframe based on duplicate variants/mutations in V_H and V_L genes"""
    df = df.copy()  # was previously modifying the original dataframe
    new_df = (
        df.groupby(["Variant", "Concentration", "Bin"])["Count"].sum().reset_index()
    )
    new_df.drop_duplicates(
        subset=["Variant", "Concentration", "Bin"], keep="first", inplace=True
    )
    new_df.rename(columns={"Count": "Rijk"}, inplace=True)
    df.drop(columns=["Barcode", "ref_barcode"], inplace=True)
    df.drop_duplicates(
        subset=["Variant", "Concentration", "Bin"], keep="first", inplace=True
    )
    for row in new_df.itertuples():
        df.loc[
            (df["Variant"] == row.Variant)
            & (df["Concentration"] == row.Concentration)
            & (df["Bin"] == row.Bin),
            "Rijk",
        ] = row.Rijk
    # df = df.sort_values(by=['Barcode', 'Variant'])
    # put together all collapsed in with combined
    return df


def add_unique_var_num(df, samples):
    """Add a unique variant number to each variant in the dataframe"""
    mut_dict = df["Variant"].value_counts().to_dict()
    # iterate through variants and assign a unique number only if it is in all samples
    accepted_variants = {
        mut: i
        for i, (mut, count) in enumerate(mut_dict.items())
        if count == len(samples) - 1
    }
    df["Variant #"] = df["Variant"].apply(
        lambda x: accepted_variants[x] if x in accepted_variants else np.nan
    )
    return df


def add_unique_barcode_num(df, samples):
    """Add a unique barcode number to each barcode in the dataframe"""
    bc_dict = df["Barcode"].value_counts().to_dict()
    # iterate through barcodes and assign a unique number only if it is in all samples
    accepted_barcodes = {
        bc: i
        for i, (bc, count) in enumerate(bc_dict.items())
        if count == len(samples) - 1
    }
    df["Barcode #"] = df["Barcode"].apply(
        lambda x: accepted_barcodes[x] if x in accepted_barcodes else np.nan
    )
    return df


def main(argv):
    start = time.time()
    params, samples = parse_args_and_read_config_match(argv)

    if not os.path.exists(params.output_dir):
        os.makedirs(params.output_dir)

    print("Reading barcodes...")
    barcodes = pd.read_csv(params.map_file_path)
    filenames, ref_freqs = matching(params, samples, barcodes)

    if len(samples) > 1:
        print("Combining barcodes...")
        to_combine = read_barcode_files(filenames, params.limit_csv_path, ref_freqs)
        combined = pd.concat(to_combine)
        # adding unique barcode number
        #         combined = add_unique_barcode_num(combined, samples)
        # Removing nM from concentration
        combined["Concentration"] = combined["Concentration"].apply(
            lambda x: re.sub(r"([0-9.]+)nM", r"\1", x) if x is not np.nan else np.nan
        )
        collapsed = output_collapsed_counts(combined)
        # adding unique variant number
        #         collapsed = add_unique_var_num(collapsed, samples)
        collapsed.dropna(subset=["Variant"], inplace=True)
        collapsed.drop("Count", axis=1, inplace=True)

        combined.rename(columns={"Count": "Rijk"}, inplace=True)
        # combined = mle.calc_er_barcodes(combined)
        # collapsed = mle.calc_er(collapsed)

        # sort by concentration and bin
        # NOTE: Could this be done at beginning? I don't think so because of output_collapsed_counts. But maybe.
        collapsed = collapsed.sort_values(
            by=["Concentration", "Bin"], ascending=[True, False]
        )
        combined = combined.sort_values(
            by=["Concentration", "Bin"], ascending=[True, False]
        )

        collapsed.to_csv(
            os.path.join(params.output_dir, "{}_collapsed.csv".format(params.gene)),
            index=False,
        )
        combined.to_csv(
            os.path.join(params.output_dir, f"{params.gene}_combined.csv"), index=False
        )
    else:
        print(
            "No need to combine barcodes. Appending reference frequencies and other data..."
        )
        file = read_barcode_files(
            [os.path.join(params.output_dir, "ref_matched.csv")],
            params.limit_csv_path,
            ref_freqs,
        )[0]
        # NOTE: should I name this reference or something else?
        file.to_csv(
            os.path.join(params.output_dir, f"{params.gene}_reference.csv"), index=False
        )
    end = time.time()
    print("Matching output written to {}".format(params.output_dir))
    print("Finished in", round((end - start), 2), "seconds.")


if __name__ == "__main__":
    main(sys.argv[1:])
    # doctest.testmod() # Comment when you actually want to use this script, uncomment to use doctests
    sys.exit(0)

# FIXME: update scan_output.csv example (just delete that last column)
# FIXME: Bring back the required files for 4A8, CC12.1