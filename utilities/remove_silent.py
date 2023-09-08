import re
import pandas as pd
from utilities.dna import is_dna_seq, translate_sequence
from utilities.mutation import Mutation, NontargetMutation
from utilities.tile import Tile
from sklearn.preprocessing import MultiLabelBinarizer


def remove_silent(mut_string):
    """Remove silent mutations from a string of mutations that have been parsed from a Mutations dataclass.
    This will be mainly useful as a utility function if you want to remove silent mutations from a string of mutations."""
    mut_list = mut_string.split("|")
    mut_list = [mut.split(";") for mut in mut_list]
    hlist = []
    for i, mut in enumerate(mut_list[0]):
        if "WT" in mut:
            hlist.append(mut)
            break
        match = re.search(
            r"(?P<gene>\w+>V[HL]:)*(?P<str>(?P<original>[A-Z])\d+(?P<mutant>[A-Z*])-[ACGT]{3})",
            mut,
        )
        original_aa = match.group("original")
        mutant_aa = match.group("mutant")
        gene = match.group("gene")
        if gene and mutant_aa:  # verify matches are happening
            hlist.append(f"{match.group('gene')}{match.group('str')}") if (
                original_aa != mutant_aa
            ) else hlist.append(f"{match.group('gene')}")
        elif mutant_aa:  # if gene string does not get matched
            hlist.append(f"{match.group('str')}") if (
                original_aa != mutant_aa
            ) else None
    # check if silent mutation is only mutation, so hlist will only have the gene match
    if len(hlist) < 2:
        hlist[0] = hlist[0] + "WT"

    llist = []
    for i, mut in enumerate(mut_list[1]):
        if "WT" in mut:
            llist.append(mut)
            break
        match = re.search(
            r"(?P<gene>\w+>V[HL]:)*(?P<str>(?P<original>[A-Z])\d+(?P<mutant>[A-Z*])-[ACGT]{3})",
            mut,
        )
        original_aa = match.group("original")
        mutant_aa = match.group("mutant")
        gene = match.group("gene")
        if gene and mutant_aa:  # verify matches are happening
            llist.append(f"{match.group('gene')}{match.group('str')}") if (
                original_aa != mutant_aa
            ) else llist.append(f"{match.group('gene')}")
        elif mutant_aa:  # if gene string does not get matched
            llist.append(f"{match.group('str')}") if (
                original_aa != mutant_aa
            ) else None
    # check if silent mutation is only mutation, so llist will only have the gene match
    if len(llist) < 2:
        llist[0] = llist[0] + "WT"

    mut_list = [hlist, llist]
    mut_string = "|".join([";".join(mut) for mut in mut_list])
    return mut_string


def convert_to_one_hot(muts, keep_mlb=False):
    """Convert all strings in a list of variant mutations to one-hot encoding.
    muts: a list of strings of mutations as parsed from a Mutations dataclass in collapse_mutations in merge_vhvl.py
    keep_mlb: if True, return the MultiLabelBinarizer object used to convert the mutations to one-hot encoding
    """
    muts = [mut.split("|") for mut in muts]
    all_muts = []
    if len(muts[0]) == 1:  # if not in VH|VL format
        for mut in muts:
            mut = mut[0]
            mut = mut.split(":")[1]
            mut = mut.split(";")
            mut = [x for x in mut if x != "WT"]
            all_muts.append(mut)
    else:
        for mut in muts:
            mut[0] = mut[0].split(":")[1]
            mut[1] = mut[1].split(":")[1]
            mut[0] = mut[0].split(";")
            mut[1] = mut[1].split(";")
            mut = mut[0] + mut[1]
            mut = [x for x in mut if x != "WT"]
            all_muts.append(mut)

    mlb = MultiLabelBinarizer(sparse_output=True)

    onehots = mlb.fit_transform(all_muts)
    if keep_mlb:
        return onehots, mlb
    else:
        return onehots


def remove_non_encoded(df, wt_csv):
    """Remove rows from a dataframe that were not encoded.
    df: a dataframe with a column named "Variant" containing the variant mutations
    wt_csv: a dataframe containing the wild type sequences for each tile, along with the encoded positions for each gene.
    """
    for gene in wt_csv["gene"].unique():
        try:
            encoded = wt_csv[wt_csv["gene"] == gene]["encoded_positions"]
        except KeyError:
            raise KeyError(
                f"Could not find column 'encoded_positions' in the WT csv. \
                Please supply encoded positions or turn off the 'encoded_positions' option in the barcode config."
            )
        # check if encoded positions are NaN: if so, do nothing
        if any(encoded.isna()):
            return df
        try:
            encoded = ";".join([encoded.iloc[i] for i in range(len(encoded))])
            encoded = set(encoded.split(";"))
        except AttributeError:
            raise AttributeError("Encoded positions must be a string.")
        df["Gene"] = [variant.split(">")[0] for variant in df["Variant"]]
        all_gene_muts = df[df["Gene"] == gene]["Variant"]
        if len(all_gene_muts) == 0:
            raise ValueError(
                f"Could not find any mutations for gene {gene} in the dataframe. \
                Check that your dataframe contains a column named 'Variant' with the variant mutations."
            )
        mlb = convert_to_one_hot(all_gene_muts, keep_mlb=True)[1]
        undesired_mutations = [
            x for x in mlb.classes_ if x.split("-")[0] not in encoded
        ]
        for mut in undesired_mutations:
            df = df.loc[~((df["Gene"] == gene) & (df["Variant"].str.contains(mut)))]
    return df


def mutations_in_seq_nosilent(tile, seq):
    """Find all non-silent mutations contained in a given sequence.

    tile: a Tile object
    seq: a DNA sequence with the same length as the tile

    Returns a tuple of Mutation/NontargetMutation objects (empty tuple
    if sequence is wild type).

    Mutations within the tile's CDS are encoded as Mutation objects
    specifying an amino acid change and codon, while mutations outside
    the tile's CDS are encoded as NontargetMutation objects specifying
    a DNA sequence change. Note that this does not currently take into
    account the Tile's list of positions.
    """
    if not is_dna_seq(seq):
        raise TypeError("seq is not a DNA sequence.")
    if len(seq) != tile.length:
        raise ValueError("seq has a different length than tile.")
    muts = []

    for i in range(tile.cds_start, tile.cds_end, 3):
        wt_codon = tile.wt_seq[i : i + 3]
        codon = seq[i : i + 3]
        if codon != wt_codon:
            pos = (i - tile.cds_start) // 3 + tile.first_aa
            wt_aa = translate_sequence(wt_codon)
            aa = translate_sequence(codon)
            if wt_aa == aa:
                continue
            muts.append(Mutation(int(pos), wt_aa, aa, codon))
    return tuple(muts)