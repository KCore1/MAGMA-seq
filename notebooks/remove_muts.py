import pandas as pd
import numpy as np
import gzip
import sys

# ASCII values of DNA characters.
bA = ord("A")
bC = ord("C")
bG = ord("G")
bT = ord("T")
bN = ord("N")


def reverse_complement(s):
    t = np.zeros_like(s)
    v = t[::-1]
    v[s == bA] = bT
    v[s == bC] = bG
    v[s == bG] = bC
    v[s == bT] = bA
    v[s == bN] = bN
    return t


def str_to_byte_array(s):
    return np.array(bytearray(s, "ascii"), dtype=np.int8)


def byte_array_to_str(a):
    return a.tobytes().decode("ascii")


def open_by_extension(path, mode):
    """Open a file using gzip.open if its name ends with '.gz', otherwise
    use open."""
    return (gzip.open if path.endswith("gz") else open)(path, mode)


def read_line(f):
    """Read a line from f with trailing whitespace stripped.

    Unlike f.readline(), this function raises an EOFError instead of returning
    an empty string at the end of the file.
    """
    line = f.readline()
    if len(line) == 0:
        raise EOFError
    return line.rstrip()


def read_seqs(f):
    """Generate FASTQ records as tuples from an open file handle."""
    while True:
        # Read the sequence ID. If there's nothing to read, then we're done.
        try:
            seq_id = read_line(f)
        except EOFError:
            return

        # If we successfully read a sequence ID, then running out of stuff to
        # read means a truncated record.
        try:
            seq = str_to_byte_array(str(read_line(f)))
            qual_id = read_line(f)
            qual = str_to_byte_array(str(read_line(f)))
        except EOFError:
            raise EOFError("EOF while reading sequence.")

        # Some simple checks of the data.
        if seq_id[0] != "@":
            raise ValueError("Sequence ID doesn't begin with '@'.")
        if qual_id[0] != "+":
            raise ValueError("Quality ID doesn't begin with '+'.")
        if len(seq) != len(qual):
            raise ValueError("Sequence and quality are different lengths.")

        yield (seq_id, seq, qual_id, qual)


### Modify these parameters
WT_SEQ = "agaaccagttctccctgaagctgagctctgtgaccgccgcagacacggctgtctattactgtgcgagacaatggaaatggttcggggaagcctggtacttcgatctctggggccgtggcaccctggtcactgtctcctcaacc".upper()
f1 = open_by_extension("../Test2/MBK-2_S2_L001_R1_001.fastq.gz", "rt")
f2 = open_by_extension("../Test2/MBK-2_S2_L001_R2_001.fastq.gz", "rt")
output_name = "merge_test_COV2.csv"
fwd_match = "TGCTAGCGTTTTAGCAGGTcagctgcagctgcaggagtcgggcccaggactggt".upper()

i = 0
ids, seqs, quals, amplens = [], [], [], []

for i, (s1, s2) in enumerate(zip(read_seqs(f1), read_seqs(f2))):
    if i % 500 == 0:
        print(i)
    id_val1, seq1, space1, qual1 = s1
    id_val2, seq2, space2, qual2 = s2
    seq_str1 = byte_array_to_str(seq1)[:248]
    rev_read = byte_array_to_str(seq2)[:36]
    if seq_str1[:54] == fwd_match:
        fwd_read = seq_str1 + WT_SEQ
        fwd_read = byte_array_to_str(reverse_complement(str_to_byte_array(fwd_read)))
        if min(qual1) > 10 + 33:
            merged = rev_read + fwd_read
            ids.append(id_val1.split(" ")[0] + str(" merged"))
            seqs.append(merged)
            quals.append(qual1)
            amplens.append(len(merged))
        else:
            continue
    i += 1

output_df = pd.DataFrame({"ID": ids, "Seq": seqs, "Qual": quals, "Amplen": amplens})
output_df.to_csv(output_name, index=False)
