"""Merge V_h and V_l genes from separate barcode-variant maps into a single map."""

import numpy as np
import pandas as pd
import ast

# from Bio.Seq import Seq
# from Bio import pairwise2
# from Bio.pairwise2 import format_alignment
from utilities.tile import Tile, mutations_in_seq
from utilities.remove_silent import mutations_in_seq_nosilent

# from pathlib import Path
import argparse


def merge_vhvl(v_h: pd.DataFrame, v_l: pd.DataFrame):
    """Merge multiple genes from separate barcode-variant maps into a single map.

    v_h, v_l: two DataFrames that contain the barcode-variant map information for V_H and V_L genes, respectively.
    These should be the output from two samples inputted into the config for a barcode-mapping run.

    Returns: A DataFrame with columns 'Barcode', 'V_H', 'V_L', 'Count', 'Freq', 'P-value', 'Gene', 'Variant' not necessarily in that order concatenated for both.
    """

    # Rename columns
    v_h = v_h.rename(
        columns={
            "Barcode": "Barcode",
            "Variant": "Variant_H",
            "Adjusted Count": "Adj_Count_H",
            "Frequency": "Frequency_H",
            "Count": "Count_H",
            "Gene": "Gene_H",
        }
    )

    v_l = v_l.rename(
        columns={
            "Barcode": "Barcode",
            "Variant": "Variant_L",
            "Adjusted Count": "Adj_Count_L",
            "Frequency": "Frequency_L",
            "Count": "Count_L",
            "Gene": "Gene_L",
        }
    )
    
    # Throw error if dataframes are empty
    if v_h.empty:
        raise ValueError("V_H DataFrame is empty.")
    if v_l.empty:
        raise ValueError("V_L DataFrame is empty.")    

    # Merge the two DataFrames
    df = pd.merge(v_h, v_l, on="Barcode", how="inner")

    df["Gene_H"] = df["Gene_H"].apply(lambda x: ast.literal_eval(x))
    df["Gene_L"] = df["Gene_L"].apply(lambda x: ast.literal_eval(x))

    # Sort out barcodes that did not get matched VH and VL to the same gene.
    failure_df = df[df.apply(lambda x: x["Gene_H"][0] != x["Gene_L"][0], axis=1)]
    df = df[df.apply(lambda x: x["Gene_H"][0] == x["Gene_L"][0], axis=1)]

    df["Variant"] = df["Variant_H"].values + "|" + df["Variant_L"].values

    return df, failure_df


# Collapse mutations and distinguish separate genes


def collapse_mutations(
    df: pd.DataFrame, wts: pd.DataFrame, args: argparse.Namespace
) -> pd.DataFrame:
    """Collapse mutations into one-hot type encoding (e.g. CR6261>VH:S31A-CODON;S32A-CODON|VL:F74S-CODON)

    df: DataFrame with columns 'Barcode', 'Variant', 'Count', 'Freq', 'P-value', 'Gene'.
    wts: dict of gene names matched with wild type sequences for VH and VL genes.
    args: barcode parameters/arguments/namespace from argparse.

    Returns: a DataFrame with columns 'Barcode', 'Variant', 'Mutations', 'Count', 'Freq', 'P-value'.
    """
    for i, item in df.iterrows():
        # Use the gene name to find the wild type sequence
        wt = wts[
            (wts["gene"] == item["Gene"][0]) & (wts["vh_or_vl"] == item["Gene"][1])
        ].iloc[0]
        tileargs = [1, 0, len(wt.seq)]
        tile = Tile(wt.seq.upper(), *tileargs)

        if args.mutations == "non-silent":
            muts = mutations_in_seq_nosilent(
                tile, item.loc[args.vh_or_vl].upper()[: len(wt.seq)]
            )
        elif args.mutations == "all":
            muts = mutations_in_seq(
                tile, item.loc[args.vh_or_vl].upper()[: len(wt.seq)]
            )
        else:
            raise ValueError(
                "Invalid mutations argument. Must be 'all', 'non-silent', or a list of mutations."
            )

        # Collapse mutations
        muts = [repr(mut) for mut in muts]
        if len(muts) == 0:
            mutations = f"{''.join(wt['gene'])}>{args.vh_or_vl}:WT"
        else:
            mutations = f"{''.join(wt['gene'])}>{args.vh_or_vl}:{';'.join(muts)}"
        df.at[i, "Variant"] = mutations

    return df


if __name__ == "__main__":
    import sys

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Parsing command line arguments for V_H and V_L barcode-variant maps."
    )
    parser.add_argument("v_h", type=str, help="Path to the barcode-variant map for V_H")
    parser.add_argument("v_l", type=str, help="Path to the barcode-variant map for V_L")
    parser.add_argument("output", type=str, help="Path to the output file")
    args = parser.parse_args()

    # Read the input files
    v_h = pd.read_csv(args.v_h)
    v_l = pd.read_csv(args.v_l)

    # Merge the two DataFrames
    df, failure_df = merge_vhvl(v_h, v_l)

    # Collapsing of mutations is already done in the haplotyping script

    # Write the output
    print(f"Writing VH and VL merged output to {args.output}")
    df.to_csv(args.output, index=False)
    failure_df.to_csv(args.output + ".failed", index=False)

    sys.exit(0)
