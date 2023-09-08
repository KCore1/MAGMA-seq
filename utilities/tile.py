import dataclasses
from dataclasses import dataclass
from typing import List

from dna import is_dna_seq, translate_sequence
from mutation import Mutation, NontargetMutation


@dataclass(frozen=True)
class Tile:
    wt_seq: str
    first_aa: int
    cds_start: int
    cds_end: int
    cds_length: int = dataclasses.field(init=False)
    positions: List[int] = None

    def __post_init__(self):
        if not is_dna_seq(self.wt_seq):
            raise TypeError("wt_seq is not a DNA sequence.")
        if not self.cds_start >= 0:
            raise ValueError("cds_start cannot be less than zero.")
        if not self.cds_end >= self.cds_start:
            raise ValueError("cds_end cannot be less than cds_start.")
        if not self.cds_end <= len(self.wt_seq):
            raise ValueError("cds_end cannot be greater than wt_seq length.")
        if (self.cds_end - self.cds_start) % 3 != 0:
            raise ValueError("cds_end - cds_start must be divisible by 3.")
        cds_seq = self.wt_seq[self.cds_start : self.cds_end]
        wt_aa = dict(enumerate(translate_sequence(cds_seq), self.first_aa))
        if self.positions is not None:
            try:
                positions = tuple(sorted(self.positions))
            except TypeError:
                raise TypeError(f"Not a valid iterable: {self.positions}")
            for p in positions:
                if p not in wt_aa:
                    raise ValueError(f"Not a valid position: {p}")
        else:
            positions = tuple(sorted(wt_aa))
        # Because this is frozen, we need to use object.__setattr__
        # instead of simple assignment.
        object.__setattr__(self, "length", len(self.wt_seq))
        object.__setattr__(self, "cds_length", self.cds_end - self.cds_start)
        object.__setattr__(self, "wt_aa", wt_aa)
        object.__setattr__(self, "positions", positions)


def mutations_in_seq(tile, seq):
    """Find all mutations contained in a given sequence.

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
            muts.append(Mutation(int(pos), wt_aa, aa, codon))

    return tuple(muts)
