import os
import gzip
from utilities.merge import (
    byte_array_to_str,
    read_seqs,
    merge_reads,
    compare_seq_ids,
    str_to_byte_array,
    reverse_complement,
)
from utilities.remove_silent import remove_non_encoded
import numpy as np
import pandas as pd
import time
from scipy.stats import binom
import warnings
from merge_vhvl import collapse_mutations

bN = ord("N")  # Taken from merge.py

mixed_bases = {
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

# ASCII value equivalent to a quality score of 0.
MIN_QUAL = 33


def open_by_extension(path, mode):
    """Open a file using gzip.open if its name ends with '.gz', otherwise
    use open."""
    return (gzip.open if path.endswith("gz") else open)(path, mode)


# For these functions, s can be string or byte array


def barcode(s, barcode_idx_start, barcode_idx_end):
    # NOTE: This does not include the last index.
    return s[barcode_idx_start:barcode_idx_end]


def cds(s, gene_start, gene_end):
    # NOTE: This does not include the last index.
    return s[gene_start:gene_end]


def is_good_barcode(read, template, verbose=False):
    """Checks the quality and template match of a given barcode.

    read: barcode segment of a merged sequencing read.
    template: template of barcode in mixed base format.
    verbose: print outputs and intermediate steps."""
    is_good = True
    errors = 0
    if len(template) > 0:
        assert len(read) == len(template)
        for i, mixed_base in enumerate(template):
            if verbose:
                print(mixed_base, mixed_bases[mixed_base], read[i])
            if read[i] not in mixed_bases[mixed_base]:
                is_good = False
                errors += 1
    return is_good, errors


def choose_variant(bc: str, d: dict, freq: dict, barcode_params, alpha=0.05):
    """Choose the variant that best fits the barcode.
    This version simply chooses the variant with the highest counts.

    bc: barcode
    d: dictionary of possible variant matches from match_barcode_and_variant for barcode bc
    freq: frequencies of possible variants, also from match_barcode_and_variant
    template: template of barcode in mixed base format
    alpha: p-value at which to accept a given barcode (NOT USED HERE)"""
    count_filter_threshold = barcode_params.count_threshold
    frequency_filter_threshold = barcode_params.frequency_threshold
    template = barcode_params.barcode
    if template != "":
        if is_good_barcode(bc, template)[0]:
            d = {key: val for key, val in d.items() if val > count_filter_threshold}
            d_adjusted = {
                i: d[i] / freq[i]
                for i in d.keys()
                if freq[i] > frequency_filter_threshold
            }
            if len(d_adjusted.keys()) == 0:
                return None, None, None
            var = max(d_adjusted, key=d_adjusted.get)
            count = round(d_adjusted[var], 1)
            alpha = d[var]
            return var, count, alpha
        else:
            return None, None, None
    else:
        d = {key: val for key, val in d.items() if val > count_filter_threshold}
        d_adjusted = {
            i: d[i] / freq[i] for i in d.keys() if freq[i] > frequency_filter_threshold
        }
        if len(d_adjusted.keys()) == 0:
            return None, None, None
        var = max(d_adjusted, key=d_adjusted.get)
        count = round(d_adjusted[var], 1)
        alpha = d[var]
        return var, count, alpha


def merge_all_reads(
    f1,
    f2,
    amplen,
    barcode_start,
    barcode_end,
    gene_start,
    gene_end,
    max_mm=None,
    min_qual=None,
    barcode_min_qual=None,
    gene_min_qual=None,
):
    """Merge fixed-length paired-end reads from two FASTQ files.

    f1: open file handle for forward read (CODING SEQUENCE)
    f2: open file handle for reverse read
    amplen: length of amplicon
    max_mm: maximum number of mismatches allowed
    min_qual: reads with any quality score lower than this are discarded
    """

    for i, (r1, r2) in enumerate(zip(read_seqs(f1), read_seqs(f2)), 1):
        seq_id1, seq1, qual_id1, qual1 = r1
        seq_id2, seq2, qual_id2, qual2 = r2

        match, prefix = compare_seq_ids(seq_id1, seq_id2)
        if not match:
            warnings.warn(f"Sequence IDs do not match: {seq_id1} {seq_id2}")
            continue  # NOTE: This is a change from the original merge_all_reads function
            # raise ValueError('Reads do not appear to match.')
        seq_id = prefix + " merged"

        s, q, n_mm = merge_reads(seq1, seq2, qual1, qual2, amplen)

        # Discard reads with length errors
        if s is None:
            warnings.warn(f"Merge failed: check amplicon length.")
            continue

        # Discard merged reads with too many mismatches.
        if max_mm is not None and n_mm > max_mm:
            warnings.warn(f"Merge failed: too many mismatches ({n_mm}).")
            continue

        # Discard merged reads with overall low quality scores.
        if min_qual and q.min() < min_qual + MIN_QUAL:
            warnings.warn(f"Merge failed: low quality score.")
            continue

        # Discard merged reads with low quality scores in barcode region.
        if (
            barcode_min_qual
            and q[barcode_start:barcode_end].min() < barcode_min_qual + MIN_QUAL
        ):
            warnings.warn(f"Merge failed: low barcode quality score.")
            continue

        # Discard merged reads with low quality scores in gene region.
        if gene_min_qual and q[gene_start:gene_end].min() < gene_min_qual + MIN_QUAL:
            warnings.warn(f"Merge failed: low gene quality score.")
            continue

        # Discard merged reads containing Ns.
        if bN in s:
            continue

        # s is encoded as a byte array. Convert it to a string before returning. NOTE: UNDID
        yield s, q, seq_id  # Quality is also a byte array here — will this matter?


def min_hamming_distance(mut, wts, print_align=False):
    distances = []
    for wt in wts:
        gene_start = wt[
            1
        ]  # NOTE: Make sure to check these if things are not working out well
        gene_end = wt[2]
        wt = wt[0].upper()  # NOTE: depends on the format of the input
        str1 = byte_array_to_str(
            reverse_complement(str_to_byte_array(cds(mut, gene_start, gene_end)))
        )
        if len(str1) != len(wt):
            distances.append(len(str1) + len(wt))  # just adding a large distance
            continue
        distances.append(sum(c1 != c2 for c1, c2 in zip(str1, wt)))
    # set max number of mutations allowed from WT sequence
    if min(distances) > 10:
        return None, None
    else:
        min_idx = distances.index(min(distances))
    return min_idx, wts[min_idx]


def match_barcode_and_variant(reads, amplens, wt_csv_path, barcode_params):
    """Collect all variants that could be associated with a given barcode.
    Used to pass to choose_variant.

    reads: the sequencing reads from a merged set of sequencing files.
    amplens: the amplicon lengths of the reads.
    wt_csv_path: path to the CSV file containing the wildtype sequences.
    barcode_params: Namespace object from argparse containing the parameters in the [Barcode] section of a config.
    """
    counts = dict()
    freq = dict()  # frequency of each variant
    genes = dict()  # gene associated with each variant
    wts = pd.read_csv(wt_csv_path)
    # at this point this will already have been verified to exist
    vh_or_vl = barcode_params.vh_or_vl
    for r, a in zip(reads, amplens):
        if vh_or_vl == "VL":  # reverse complement of gene
            # possible_wts = [wt for wt in wts.itertuples(
            #     index=False) if (wt.amplen == a and wt.vh_or_vl == 'VL')]
            possible_wts = [
                wt for wt in wts.itertuples(index=False) if (wt.vh_or_vl == "VL")
            ]
            possible_VLs = [(wt.seq, wt.gene_start, wt.gene_end) for wt in possible_wts]
            if len(possible_wts) == 0:
                raise ValueError(
                    "No wildtype sequences found for amplicon length " + str(a)
                )
            wt_idx = min_hamming_distance(r, possible_VLs)[0]
            if wt_idx is None:
                continue
            else:
                wt = possible_wts[wt_idx]
                start = wt.barcode_start
                end = wt.barcode_end
                gene_start = wt.gene_start
                gene_end = wt.gene_end
                var = byte_array_to_str(
                    reverse_complement(str_to_byte_array(cds(r, gene_start, gene_end)))
                )
                bc = barcode(r, start, end)
        else:
            # possible_wts = [wt for wt in wts.itertuples(
            #     index=False) if (wt.amplen == a and wt.vh_or_vl == 'VH')]
            possible_wts = [
                wt for wt in wts.itertuples(index=False) if (wt.vh_or_vl == "VH")
            ]
            possible_VHs = [(wt.seq, wt.gene_start, wt.gene_end) for wt in possible_wts]
            if len(possible_wts) == 0:
                raise ValueError(
                    "No wildtype sequences found for amplicon length " + str(a)
                )
            wt_idx = min_hamming_distance(r, possible_VHs)[0]
            if wt_idx is None:
                continue
            else:
                wt = possible_wts[wt_idx]
                start = wt.barcode_start
                end = wt.barcode_end
                gene_start = wt.gene_start
                gene_end = wt.gene_end
                var = byte_array_to_str(
                    reverse_complement(str_to_byte_array(cds(r, gene_start, gene_end)))
                )
                bc = byte_array_to_str(
                    reverse_complement(str_to_byte_array(barcode(r, start, end)))
                )
        counts.setdefault(bc, dict())
        counts[bc].setdefault(var, 0)  # This will be a dict
        freq.setdefault(var, 0)
        counts[bc][var] += 1
        freq[var] += 1
        genes[var] = (wt.gene, wt.vh_or_vl)

    total_size = sum(freq.values())
    freq = {k: v / total_size for k, v in freq.items()}  # normalize frequencies

    return counts, freq, genes


def barcode_to_variant_map(reads, amplens, wt_csv_path, barcode_params):
    """Maps barcodes to gene variant, using the choose_variant, match_barcode_and_variant functions.
    Requires merged and filtered reads (in str format) as input. Returns three dicts (one of successfully
    chosen, one of failed, one of all the read variants) of barcodes to variants data, including counts and p-values.
    """

    counts, freq, genes = match_barcode_and_variant(
        reads, amplens, wt_csv_path, barcode_params
    )

    # bc_var_pairs_success = [(bc, *choose_variant(bc, d, freq, barcode_params.barcode), freq[choose_variant(bc, d, freq, barcode_params.barcode)[0]])
    # for bc, d in counts.items() if choose_variant(bc, d, freq, barcode_params.barcode)[0] is not None] # Fancy list comprehension but I think it's more readable in the loop
    bc_var_pairs_success = []  # Comment if above is used
    bc_var_pairs = []
    bc_var_pairs_failure = []

    for bc, d in counts.items():
        bc_var_pairs.extend(
            [
                (
                    bc,
                    mut,
                    count,
                    1,
                    freq[mut],
                    genes[mut],
                )  # 1 is hack for now to get this to work, talk to Brian about this
                for mut, count in d.items()
            ]
        )
        chosen_variant, count, p_val = choose_variant(bc, d, freq, barcode_params)
        if chosen_variant is not None:
            bc_var_pairs_success.append(
                (
                    bc,
                    chosen_variant,
                    count,
                    p_val,
                    freq[chosen_variant],
                    genes[chosen_variant],
                )
            )
        else:
            bc_var_pairs_failure.extend(
                [
                    (bc, mut, count_d, 1, freq[mut], genes[mut])
                    for mut, count_d in d.items()
                ]
            )  # also has 1

    bc_map_success = pd.DataFrame(
        bc_var_pairs_success,
        columns=[
            "Barcode",
            barcode_params.vh_or_vl,
            "Adjusted Count",
            "Count",
            "Frequency",
            "Gene",
        ],
    )
    bc_map_failure = pd.DataFrame(
        bc_var_pairs_failure,
        columns=[
            "Barcode",
            barcode_params.vh_or_vl,
            "Adjusted Count",
            "Count",
            "Frequency",
            "Gene",
        ],
    )
    bc_map = pd.DataFrame(
        bc_var_pairs,
        columns=[
            "Barcode",
            barcode_params.vh_or_vl,
            "Adjusted Count",
            "Count",
            "Frequency",
            "Gene",
        ],
    )

    return bc_map_success, bc_map_failure, bc_map


def write_barcode_map(filename, barcode_map, debug=False):
    if debug:
        start = time.time()
    barcode_map.to_csv(filename, index=False)
    if debug:
        end = time.time()
        print(f"Writing barcode map to {filename} took {end-start} seconds.")


def write_merged_reads(
    filepath,
    path1,
    path2,
    wt_csv_path,
    vh_or_vl,
    max_mm,
    min_qual,
    barcode_min_qual,
    gene_min_qual,
    debug: bool = False,
):
    """Writes merged reads to a file. Requires open file handles for f1 and f2. Returns the dataframe of merged reads.

    filepath: Path to write merged reads to.
    f1, f2: open file handles to sequences in fastq or fastq.gz format.
    wt_csv_path: path to csv containing information about each gene and its wild type sequence.
    max_mm: max number of mismatches allowed in overlap region.
    barcode_min_qual: minimum quality allowed in barcode region.
    barcode_start: starting index of the barcode.
    barcode_end: ending index of the barcode (not inclusive).
    debug: turn on or off timing."""
    ids, seqs, quals, amplens = (
        [],
        [],
        [],
        [],
    )  # FUTURE: Pre allocate np arrays for vectorized operations when adding them all to df

    wts = pd.read_csv(wt_csv_path)
    amplens_to_merge_at = set(wts[wts["vh_or_vl"] == vh_or_vl]["amplen"])
    # make sure barcode regions are of high quality
    barcode_start = min(wts["barcode_start"])
    barcode_end = max(wts["barcode_end"])
    gene_start = min(wts["gene_start"])
    gene_end = max(wts["gene_end"])
    for amplen in amplens_to_merge_at:
        f1 = open_by_extension(path1, "rt")
        f2 = open_by_extension(path2, "rt")
        for seq, qual, id in merge_all_reads(
            f1,
            f2,
            amplen,
            barcode_start,
            barcode_end,
            gene_start,
            gene_end,
            max_mm,
            min_qual,
            barcode_min_qual,
            gene_min_qual,
        ):
            start = time.time() if debug else None
            seq_str = byte_array_to_str(seq)
            qual_str = byte_array_to_str(qual)

            ids.append(id)
            seqs.append(seq_str)
            quals.append(qual_str)
            amplens.append(amplen)

            end = time.time() if debug else None

            print(f"merged read {id} in {end - start} seconds.") if debug else None
        f1.close()
        f2.close()
        print("merged all reads at amplicon length", amplen) if debug else None

    # Concat is quadratic, so just creating a df from these lists is much faster
    start = time.time()
    df = pd.DataFrame({"ID": ids, "Seq": seqs, "Qual": quals, "Amplen": amplens})

    if os.path.exists(filepath):
        df.to_csv(filepath, index=False)
    else:
        os.makedirs(os.path.dirname(filepath))
        df.to_csv(filepath, index=False)

    end = time.time()
    print(
        f"saved merged reads to {filepath} in {end - start} seconds."
    ) if debug else None
    return df


# NOTE: Passing in params seems kinda clunky, not sure what to do otherwise though
def get_barcode_map(path1, path2, params, barcode_params, sample_name):
    pathname = params.output_dir
    wt_csv_path = barcode_params.wt_csv
    if os.path.exists(f"{pathname}/{sample_name}/merged.csv"):
        # Reading from csv for faster runtime
        print(f"Reading {sample_name} merged reads from csv...")
        df = pd.read_csv(f"{pathname}/{sample_name}/merged.csv")
    else:
        print(f"Merging reads for sample {sample_name}...")
        try:
            _ = barcode_params.vh_or_vl  # Try to access variable
        except ValueError:
            # FUTURE: Add support for other genes?
            raise ValueError("vh_or_vl parameter must be VH or VL")

        if (
            barcode_params.vh_or_vl == "VH"
        ):  # VH needs have flipped forward and reverse, then also needs to be reverse complemented
            df = write_merged_reads(
                f"{pathname}/{sample_name}/merged.csv",
                path2,
                path1,
                wt_csv_path,
                barcode_params.vh_or_vl,
                params.max_mismatches,
                params.min_quality,
                barcode_params.barcode_min_quality,
                params.gene_min_quality,
            )  # create merged file
        # VL needs to be just reverse complemented on the gene (done in match_barcode_and_variant)
        else:
            df = write_merged_reads(
                f"{pathname}/{sample_name}/merged.csv",
                path1,
                path2,
                wt_csv_path,
                barcode_params.vh_or_vl,
                params.max_mismatches,
                params.min_quality,
                barcode_params.barcode_min_quality,
                params.gene_min_quality,
            )  # create merged file
    num_lines = sum(
        1 for _ in open_by_extension(path1, "rt")
    )  # sum number of lines in the file
    num_reads_in_file = num_lines // 4  # FASTQ format has 4 lines per read
    print("Number of reads in file:", num_reads_in_file)
    print("Number of merged reads:", len(df))

    reads = df.loc[:, "Seq"]
    amplens = df.loc[:, "Amplen"]
    success_map, failure_map, bc_map = barcode_to_variant_map(
        reads, amplens, wt_csv_path, barcode_params
    )

    # extract mutations using function from from merge_vhvl.py
    wts_csv = pd.read_csv(barcode_params.wt_csv)
    print(f"Extracting mutations for sample {sample_name}...")

    # NOTE: Using version of mutations_from_seq in collapse_mutations that does not include silent mutatations
    success_map, failure_map, bc_map = (
        collapse_mutations(success_map, wts_csv, barcode_params),
        collapse_mutations(failure_map, wts_csv, barcode_params),
        collapse_mutations(bc_map, wts_csv, barcode_params),
    )

    if barcode_params.encoded_positions:
        success_map, failure_map, bc_map = (
            remove_non_encoded(success_map, wts_csv).reset_index(drop=True),
            remove_non_encoded(failure_map, wts_csv).reset_index(drop=True),
            remove_non_encoded(bc_map, wts_csv).reset_index(drop=True),
        )

    success_map.sort_values(
        by=["Barcode", "Count"], ascending=[True, False], inplace=True
    )
    failure_map.sort_values(
        by=["Barcode", "Count"], ascending=[True, False], inplace=True
    )
    bc_map.sort_values(by=["Barcode", "Count"], ascending=[True, False], inplace=True)

    print(f"Writing barcode maps for sample {sample_name}...")
    write_barcode_map(f"{pathname}/{sample_name}/map_failure.csv", failure_map)
    print(f"Done with {sample_name} failure map.")
    write_barcode_map(f"{pathname}/{sample_name}/map_success.csv", success_map)
    print(f"Done with {sample_name} success map.")
    write_barcode_map(f"{pathname}/{sample_name}/map.csv", bc_map)
    print(f"Done with {sample_name} map.")
