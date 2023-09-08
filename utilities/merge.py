import numpy as np

# ASCII values of DNA characters.
bA = ord("A")
bC = ord("C")
bG = ord("G")
bT = ord("T")
bN = ord("N")

# ASCII value equivalent to a quality score of 0.
MIN_QUAL = 33

# These are written like this for speed because this is the most
# expensive part of the code. Writing this in C would speed it up
# substantially, but would also make distribution more difficult.


# These don't do any checking that only valid bases are present, but
# it will become obvious quickly in the rest of the code if files
# contain invalid DNA sequences.
def str_to_byte_array(s):
    return np.array(bytearray(s, "ascii"), dtype=np.int8)


def byte_array_to_str(a):
    return a.tobytes().decode("ascii")


def reverse_complement(s):
    t = np.zeros_like(s)
    v = t[::-1]
    v[s == bA] = bT
    v[s == bC] = bG
    v[s == bG] = bC
    v[s == bT] = bA
    v[s == bN] = bN
    return t


def read_line(f):
    """Read a line from f with trailing whitespace stripped.

    Unlike f.readline(), this function raises an EOFError instead of returning
    an empty string at the end of the file.
    """
    line = f.readline()
    if len(line) == 0:
        raise EOFError
    return line.rstrip()


def merge_reads(
    s1, s2, q1, q2, amplen
):
    """Merge paired end reads of an amplicon and return sequence,
    quality, and number of mismatches.

    s1 -- Sequence of first read (ndarray of int).
    s2 -- Sequence of second read (ndarray of int).
    q1 -- Quality of first read (ndarray of int).
    q2 -- Quality of second read (ndarray of int).
    amplen -- Length of amplicon (int).

    """
    # If the amplicon is of length L and the reads are lengths l1, l2 then:
    # - read 1 from 0 to L-l2-1 inclusive doesn't overlap
    # - read 1 from L-l2 to l1-1 inclusive overlaps with read 2
    # - read 2 from 0 to l1+l2-L-1 inclusive overlaps with read 1
    # - read 2 from l1+l2-L to its end doesn't overlap

    # A picture for clarity:
    # s1 coords: 0                                      l1-1
    #            |                                      |
    #            ----------------------------------------
    #                                          ------------------------------
    #                                          |        |                   |
    # s1 coords:                               L-l2     |                   L-1
    # s2 coords:                               0        l1+l2-L-1

    # Reverse complement read 2 and reverse its quality scores.
    s2 = reverse_complement(s2)
    q2 = q2[::-1]

    # This is where we'll put the merged sequence and quality score.
    s = np.zeros(amplen, dtype=np.int8)
    q = np.zeros(amplen, dtype=np.int8)

    # If the reads overlap correctly, then s1[offset+i] == s2[i], assuming s2 is
    # the reverse complement of the reverse read.
    offset = amplen - len(
        s2
    )  # Expecting read lengths to be shorter than amplicon length.

    # Fill in the parts of the merged sequence where the reads don't overlap.

    if (
        offset < 0
    ):  # Make sure the lengths are correct if reverse read is longer than amplicon.
        s2 = s2[-offset:]
        q2 = q2[-offset:]
        offset = amplen - len(s2)
    if len(s2) + len(s1) < amplen:  # Make sure reads overlap
        return None, None, None
    s[:offset] = s1[:offset]
    q[:offset] = q1[:offset]
    s[len(s1) :] = s2[len(s1) + len(s2) - amplen :]
    q[len(s1) :] = q2[len(s1) + len(s2) - amplen :]

    # Create a set of views into the overlapping region. We can directly compare
    # vs1[i] to vs2[i] and use that to fill in vs[i] with all indexing taken
    # care of.
    vs1 = s1[offset:]
    vq1 = q1[offset:]
    vs2 = s2[: len(vs1)]
    vq2 = q2[: len(vs1)]
    vs = s[offset : len(s1)]
    vq = q[offset : len(s1)]

    # Quality score of matching bases is the larger of the two quality
    # scores (this is a somewhat conservative low estimate). Quality
    # score of mismatched bases is the difference of the two quality
    # scores. If the mismatched bases have equal quality scores, the
    # base is written as an N with the minimum possible quality.

    # Positions where the reads agree.
    if len(vs1) == len(
        vs2
    ):  # Make sure the lengths are correct again so that the indices can be compared
        ieq = vs1 == vs2
    else:
        return None, None, None
    vs[ieq] = vs1[ieq]
    vq[ieq] = np.maximum(vq1[ieq], vq2[ieq])

    # Positions where the reads disagree.
    ineq = vs1 != vs2
    mismatches = ineq.sum()

    # Positions where the reads disagree and read 1 has the higher quality.
    ir1 = np.logical_and(ineq, vq1 > vq2)
    vs[ir1] = vs1[ir1]
    vq[ir1] = MIN_QUAL + vq1[ir1] - vq2[ir1]

    # Positions where the reads disagree and read 2 has the higher quality.
    ir2 = np.logical_and(ineq, vq2 > vq1)
    vs[ir2] = vs2[ir2]
    vq[ir2] = MIN_QUAL + vq2[ir2] - vq1[ir2]

    # Positions where the reads disagree and they have equal qualities.
    irn = np.logical_and(ineq, vq1 == vq2)
    vs[irn] = bN
    vq[irn] = MIN_QUAL

    return s, q, mismatches


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
            seq = str_to_byte_array(read_line(f))
            qual_id = read_line(f)
            qual = str_to_byte_array(read_line(f))
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


def compare_seq_ids(id1, id2):
    """Compare two sequence IDs and return their common prefix if they
    match.

    Splitting an Illumina sequence identifier gives a prefix string
    which should match between paired end reads and a read-specific
    string. If the prefix matches, this returns (True, prefix). IF the
    prefix does not match, this returns (False, None).

    This function does not worry about if the IDs start with '@'.
    """
    split1 = id1.split(" ")
    split2 = id2.split(" ")
    if len(split1) != 2:
        raise ValueError("id1 does not contain exactly one space.")
    if len(split2) != 2:
        raise ValueError("id2 does not contain exactly one space.")
    l1, r1 = split1
    l2, r2 = split2
    match = l1 == l2
    prefix = l1 if match else None
    return match, prefix


def merge_all_reads(f1, f2, amplen, max_mm=None, min_qual=None):
    """Merge fixed-length paired-end reads from two FASTQ files.

    f1: open file handle for forward read
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
            raise ValueError("Reads do not appear to match.")
        seq_id = prefix + " merged"

        s, q, n_mm = merge_reads(seq1, seq2, qual1, qual2, amplen)

        # Discard reads with length errors
        if s is None:
            continue

        # Discard merged reads with too many mismatches.
        if max_mm is not None and n_mm > max_mm:
            continue

        # Discard merged reads with low quality scores.
        if min_qual is not None and q.min() < min_qual + MIN_QUAL:
            continue

        # Discard merged reads containing Ns.
        if bN in s:
            continue

        # s is encoded as a byte array. Convert it to a string before returning.
        yield byte_array_to_str(s)
