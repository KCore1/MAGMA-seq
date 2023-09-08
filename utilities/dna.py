import math

_BASES = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
}

_AA_CODONS = {
    "A": ("GCT", "GCC", "GCA", "GCG"),
    "R": ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
    "N": ("AAT", "AAC"),
    "D": ("GAT", "GAC"),
    "C": ("TGT", "TGC"),
    "Q": ("CAA", "CAG"),
    "E": ("GAA", "GAG"),
    "G": ("GGT", "GGC", "GGA", "GGG"),
    "H": ("CAT", "CAC"),
    "I": ("ATT", "ATC", "ATA"),
    "L": ("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
    "K": ("AAA", "AAG"),
    "M": ("ATG",),
    "F": ("TTT", "TTC"),
    "P": ("CCT", "CCC", "CCA", "CCG"),
    "S": ("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
    "T": ("ACT", "ACC", "ACA", "ACG"),
    "W": ("TGG",),
    "Y": ("TAT", "TAC"),
    "V": ("GTT", "GTC", "GTA", "GTG"),
    "*": ("TAA", "TAG", "TGA"),
}

_TRANSLATE = {}
for aa, codons in _AA_CODONS.items():
    for codon in codons:
        _TRANSLATE[codon] = aa

AMINO_ACIDS = list([aa for aa in _AA_CODONS.keys() if aa != "*"])
AMINO_ACIDS_PLUS_STOP = list([aa for aa in _AA_CODONS.keys()])


def is_dna_base(s):
    try:
        return s in _BASES
    except TypeError:
        return False


def is_dna_seq(s):
    try:
        return all([b in _BASES for b in s])
    except TypeError:
        return False


def is_codon(s):
    try:
        return len(s) == 3 and is_dna_seq(s)
    except TypeError:
        return False


def reverse_complement(s):
    return "".join(reversed([_BASES[b] for b in s]))


def is_aa(aa):
    try:
        return aa in _AA_CODONS
    except TypeError:
        return False


def is_aa_seq(s):
    try:
        return all([aa in _AA_CODONS for aa in s])
    except TypeError:
        return False


def aa_codons(aa):
    if not is_aa(aa):
        raise ValueError("Not an amino acid.")
    return _AA_CODONS[aa]


def translate_sequence(s):
    if not is_dna_seq(s):
        raise ValueError("Not a DNA sequence.")
    if len(s) % 3 != 0:
        raise ValueError("Length not a multiple of three.")
    return "".join([_TRANSLATE[s[i : i + 3]] for i in range(0, len(s), 3)])
