import dataclasses
from dataclasses import dataclass

from utilities.dna import is_aa, is_codon, is_dna_base, translate_sequence


@dataclass(frozen=True)
class Mutation:
    """An amino acid mutation in a coding sequence of interest."""

    pos: int
    wt_aa: str
    aa: str
    codon: str

    def __post_init__(self):
        if not isinstance(self.pos, int):
            raise TypeError("pos must be an integer.")
        if not self.pos > 0:
            raise ValueError("pos must be positive.")
        if not is_aa(self.wt_aa):
            raise TypeError("wt_aa must be a valid amino acid or stop.")
        if not is_aa(self.aa):
            raise TypeError("aa must be a valid amino acid or stop.")
        if not is_codon(self.codon):
            raise TypeError("codon must be a valid codon.")
        if not translate_sequence(self.codon) == self.aa:
            raise ValueError("codon does not match aa.")

    def __repr__(self):
        return "%s%i%s-%s" % (self.wt_aa, self.pos, self.aa, self.codon)


@dataclass(frozen=True)
class NontargetMutation:
    """A DNA mutation outside a coding sequence of interest."""

    pos: int
    wt_base: str
    base: str

    def __post_init__(self):
        if not isinstance(self.pos, int):
            raise TypeError("pos must be an integer.")
        if not self.pos >= 0:
            raise ValueError("pos cannot be negative.")
        if not is_dna_base(self.wt_base):
            raise TypeError("wt_base must be a DNA base.")
        if not is_dna_base(self.base):
            raise TypeError("base must be a DNA base.")
        if self.wt_base == self.base:
            raise ValueError("wt_base and base cannot be equal.")


@dataclass(frozen=True)
class _WildType(object):
    def __repr__(self):
        return "WT"


WildType = _WildType()
"""Value to represent wild type sequences."""


# This is necessary because the use of multiprocessing to process
# FASTQ files in parallel means that this file will be imported once
# per process, leading to multiple WildType objects floating around.
def is_wt(m):
    return isinstance(m, _WildType)


@dataclass(frozen=True)
class AminoAcidMutation:
    """An abstract amino acid mutation without codon information."""

    pos: int
    wt_aa: str
    aa: str

    def __post_init__(self):
        if not isinstance(self.pos, int):
            raise TypeError("pos must be an integer.")
        if not self.pos > 0:
            raise ValueError("pos must be positive.")
        if not is_aa(self.wt_aa):
            raise TypeError("wt_aa must be a valid amino acid or stop.")
        if not is_aa(self.aa):
            raise TypeError("aa must be a valid amino acid or stop.")

    def __repr__(self):
        return "%s%i%s" % (self.wt_aa, self.pos, self.aa)

    @staticmethod
    def from_mutation(mut):
        if mut.aa == mut.wt_aa:
            return WildType
        return AminoAcidMutation(pos=mut.pos, wt_aa=mut.wt_aa, aa=mut.aa)

    @staticmethod
    def sort_key(x):
        return (-1, "A", "A") if is_wt(x) else (x.pos, x.wt_aa, x.aa)
