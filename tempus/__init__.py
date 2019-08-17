"""
Model classes and constants.
"""
from dataclasses import dataclass, field
from typing import Optional

from enum import Enum, unique, IntEnum


@dataclass(frozen=True)
class SimpleVariant:
    """
    Trivial representation of a single alt allele variant.
    """
    contig: str
    pos: int
    ref: str
    alt: str
    alt_index: Optional[int] = field(default=None)


class Assembly(Enum):
    """
    Reference assembly identifiers.
    """
    GRCH37: 'Assembly' = 'GRCh37'


@unique
class Impact(IntEnum):
    """
    Represents degrees of severity.
    """
    LOW: 'Impact' = 0
    MODERATE: 'Impact' = 1
    HIGH: 'Impact' = 2


@unique
class SequenceAlteration(Enum):
    """
    A sequence_alteration is a sequence_feature whose extent is the deviation from another sequence [SO:0001059].
    Annotations roughly adhere to the http://www.sequenceontology.org/ nomenclature.
    """
    DELETION: 'SequenceAlteration' = ('SO:0000159', 'deletion')
    DELINS: 'SequenceAlteration' = ('SO:1000032', 'delins')
    DUPLICATION: 'SequenceAlteration' = ('SO:1000035', 'duplication')
    INSERTION: 'SequenceAlteration' = ('SO:0000667', 'insertion')
    INVERSION: 'SequenceAlteration' = ('SO:1000036', 'inversion')
    STRUCTURAL_ALTERATION: 'SequenceAlteration' = ('SO:0001785', 'structural_alteration')
    SUBSTITUTION: 'SequenceAlteration' = ('SO:1000002', 'substitution')

    def __init__(self, so: str, slug: str):
        self.so = so
        self.slug = slug

    def __str__(self):
        return self.slug


@unique
class FeatureVariant(Enum):
    """
    A sequence variant that falls entirely or partially within a genomic feature [SO:0001878].
    Annotations roughly adhere to the http://www.sequenceontology.org/ nomenclature.
    """
    FRAMESHIFT: 'FeatureVariant' = ('SO:0001589', 'frameshift_variant', Impact.HIGH)  # coding, tx, protein alt
    TERMINATOR: 'FeatureVariant' = ('SO:0001590', 'terminator_codon_variant', Impact.HIGH)  # coding, tx, protein alt
    SYNONYMOUS: 'FeatureVariant' = ('SO:0001819', 'synonymous_variant', Impact.LOW)  # coding, tx
    MISSENSE: 'FeatureVariant' = ('SO:0001583', 'missense_variant', Impact.MODERATE)  # (non)coding, tx, protein alt
    INTRONIC: 'FeatureVariant' = ('SO:0001627', 'intron_variant', Impact.MODERATE)  # (non)coding, tx
    EXONIC: 'FeatureVariant' = ('SO:0001791', 'exon_variant', Impact.MODERATE)  # (non)coding, tx (could be UTRs)
    INTERGENIC: 'FeatureVariant' = ('SO:0001628', 'intergenic_variant', Impact.MODERATE)

    def __init__(self, so: str, slug: str, impact: Impact):
        self.so = so
        self.slug = slug
        self.impact = impact

    def __str__(self):
        return self.slug
