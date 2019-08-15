"""
Convenience functions for working with `biocommons/hgvs` package.
"""
import logging
from dataclasses import dataclass, field
from typing import Dict, Optional

from bioutils.sequences import reverse_complement
from hgvs.assemblymapper import AssemblyMapper
from hgvs.dataproviders import uta
from hgvs.edit import NARefAlt, Dup, Inv
from hgvs.exceptions import HGVSError
from hgvs.location import Interval, SimplePosition
from hgvs.normalizer import Normalizer
from hgvs.posedit import PosEdit
from hgvs.sequencevariant import SequenceVariant

from tempus import Assembly, SimpleVariant, SequenceAlteration, FeatureVariant

#
# hgvs utility code
#

ASSEMBLY__HGVS_ASSEMBLY_NAME: Dict[Assembly, str] = {
    Assembly.GRCH37: 'GRCh37'
}


@dataclass(frozen=True)
class HgvsMachinery:
    assembly_mapper: AssemblyMapper
    normalizer_3p: Normalizer
    normalizer_5p: Normalizer
    accession__contig: Dict[str, str] = field(repr=False)
    contig__accession: Dict[str, str] = field(repr=False)

    @classmethod
    def from_assembly(cls, assembly: Assembly, alt_aln_method: str = 'splign') -> 'HgvsMachinery':
        """
        Initializes `biocommons/hgvs` machinery with Universal Transcript Archive (UTA) data provider.

        :param alt_aln_method: alignment method (default 'splign').
        :param assembly: determines the assembly build.
        """
        data_provider = uta.connect()
        assembly_name = ASSEMBLY__HGVS_ASSEMBLY_NAME[assembly]
        accession__contig = data_provider.get_assembly_map(assembly_name=assembly_name)
        return cls(
            assembly_mapper=AssemblyMapper(
                hdp=data_provider, assembly_name=assembly_name, alt_aln_method=alt_aln_method, add_gene_symbol=False,
                prevalidation_level='NONE', normalize=True, replace_reference=False),
            normalizer_3p=Normalizer(hdp=data_provider, shuffle_direction=3),
            normalizer_5p=Normalizer(hdp=data_provider, shuffle_direction=5),
            accession__contig=accession__contig,
            contig__accession={contig: acc for acc, contig in accession__contig.items()})

    def hgvs_from_simple_variant(self, variant: SimpleVariant, vtype: str) -> SequenceVariant:
        return SequenceVariant(
            ac=self.contig__accession[variant.contig],
            type=vtype,
            posedit=PosEdit(
                pos=Interval(start=SimplePosition(variant.pos), end=SimplePosition(variant.pos + len(variant.ref) - 1)),
                edit=NARefAlt(ref=variant.ref, alt=variant.alt)))

    def simple_variant_from_hgvs(self, variant: SequenceVariant) -> SimpleVariant:
        edit = variant.posedit.edit
        if isinstance(edit, Dup):
            alt = edit.ref_s * 2
        elif isinstance(edit, Inv):
            alt = reverse_complement(edit.ref_s)
        else:
            alt = edit.alt
        return SimpleVariant(
            contig=self.accession__contig[variant.ac], pos=variant.posedit.pos.start.base,
            ref=variant.posedit.edit.ref, alt=alt)

    def gene_from_transcripts(self, *tx_accessions: str) -> Optional[str]:
        candidates = {self.assembly_mapper.hdp.get_tx_identity_info(tx_ac).get("hgnc", None) for tx_ac in tx_accessions}
        gene = next(iter(candidates), None)
        if len(candidates) > 1:
            logging.warning(f'picking gene {gene} from many candidates: {candidates}')
        return gene


def variant_edit_type(variant: SequenceVariant) -> Optional[str]:
    try:
        return variant.posedit.edit.type
    except AttributeError:
        return None


def is_transcript_coding(transcript: str) -> bool:
    return transcript.startswith('NM') or transcript.startswith('XM')


def is_variant_intronic(variant_c: SequenceVariant) -> bool:
    pos: Interval = variant_c.posedit.pos
    return pos.start.is_intronic and pos.end.is_intronic


#
# annotation related code
#

_HGVS_G_EDIT_TYPE__SEQUENCE_ALTERATION: Dict[str, SequenceAlteration] = {
    'del': SequenceAlteration.DELETION,
    'delins': SequenceAlteration.DELINS,
    'dup': SequenceAlteration.DUPLICATION,
    'ins': SequenceAlteration.INSERTION,
    'inv': SequenceAlteration.INSERTION,
    'sub': SequenceAlteration.SUBSTITUTION,
}

_HGVS_P_EDIT_TYPE__FEATURE_VARIANT: Dict[str, FeatureVariant] = {
    'identity': FeatureVariant.SYNONYMOUS,
    'ext': FeatureVariant.TERMINATOR,
}


@dataclass(frozen=True)
class HgvsTranscriptAnnotation:
    hgvs_c: Optional[SequenceVariant]
    hgvs_p: Optional[SequenceVariant]
    feature_variant: Optional[FeatureVariant]


@dataclass(frozen=True)
class HgvsVariantAnnotation(HgvsTranscriptAnnotation):
    hgvs_g: SequenceVariant
    sequence_alteration: SequenceAlteration
    gene: Optional[str]


def sequence_alteration_from_hgvs_g(hgvs_g: SequenceVariant) -> SequenceAlteration:
    return _HGVS_G_EDIT_TYPE__SEQUENCE_ALTERATION.get(
        variant_edit_type(hgvs_g), SequenceAlteration.STRUCTURAL_ALTERATION)


def feature_variant_from_hgvs_c(hgvs_c: SequenceVariant) -> FeatureVariant:
    return FeatureVariant.INTRONIC if is_variant_intronic(hgvs_c) else FeatureVariant.EXONIC


def feature_variant_from_hgvs_p(hgvs_p: SequenceVariant, hgvs_c: SequenceVariant) -> FeatureVariant:
    return _HGVS_P_EDIT_TYPE__FEATURE_VARIANT.get(
        variant_edit_type(hgvs_p),
        FeatureVariant.MISSENSE if hgvs_c.posedit.length_change() % 3 == 0 else FeatureVariant.FRAMESHIFT)


def annotate_variant_hgvs(hgvs_machinery: HgvsMachinery, variant: SimpleVariant) -> HgvsVariantAnnotation:
    def annotate_transcript(tx_ac: str) -> HgvsTranscriptAnnotation:
        hgvs_c: Optional[SequenceVariant] = None
        hgvs_p: Optional[SequenceVariant] = None
        feature_variant: Optional[FeatureVariant] = None
        try:
            hgvs_c = hgvs_machinery.assembly_mapper.g_to_c(hgvs_g, tx_ac)
        except HGVSError:
            pass
        else:
            feature_variant = feature_variant_from_hgvs_c(hgvs_c)
            try:
                hgvs_p = hgvs_machinery.assembly_mapper.c_to_p(hgvs_c)
            except HGVSError:
                pass
            else:
                feature_variant = feature_variant_from_hgvs_p(hgvs_p, hgvs_c)
        return HgvsTranscriptAnnotation(hgvs_c, hgvs_p, feature_variant)

    hgvs_g = hgvs_machinery.normalizer_3p.normalize(hgvs_machinery.hgvs_from_simple_variant(variant, vtype='g'))

    # transcripts and proteins
    txs_all = hgvs_machinery.assembly_mapper.relevant_transcripts(hgvs_g)

    # determine gene
    gene: Optional[str] = hgvs_machinery.gene_from_transcripts(*txs_all)

    # annotate coding transcripts (hgvs_g -> hgvs_c -> hgvs_p)
    tx_anns = tuple(annotate_transcript(tx_ac) for tx_ac in txs_all if is_transcript_coding(tx_ac))

    # pick tx annotation of max impact or pick the first one (most deleterious)
    tx_ann: Optional[HgvsTranscriptAnnotation] = max(
        tx_anns, key=lambda ann: ann.feature_variant.impact if ann.feature_variant else -1) if tx_anns else None

    feature_variant_default: FeatureVariant = FeatureVariant.INTERGENIC if gene else FeatureVariant.INTRONIC

    return HgvsVariantAnnotation(
        hgvs_c=getattr(tx_ann, 'hgvs_c', None),
        hgvs_p=getattr(tx_ann, 'hgvs_p', None),
        feature_variant=getattr(tx_ann, 'feature_variant', feature_variant_default),
        hgvs_g=hgvs_g,
        sequence_alteration=sequence_alteration_from_hgvs_g(hgvs_g),
        gene=gene)
