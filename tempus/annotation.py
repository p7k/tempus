"""
Top-level annotation functions and serialization utils.
"""
import csv
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, Dict, Any, TextIO

from flatten_dict import flatten
from hgvs.sequencevariant import SequenceVariant
from vcf import Reader
from vcf.model import _Record

from tempus import FeatureVariant, SequenceAlteration, SimpleVariant
from tempus.exac import ExacVariantAnnotation
from tempus.hgvs import HgvsVariantAnnotation, HgvsMachinery
from tempus.vcf import simple_variants_from_record, TEMPUS_REFERENCE__ASSEMBLY, VcfVariantAnnotation


@dataclass(frozen=True)
class VariantAnnotation:
    """
    Locus-level annotation encapsulating the following:

        var :: chrom-pos-ref-alt form of variant
        vcf :: annotations extracted directly from the vcf itself
        hgvs :: annotations acquired via UTA (RefSeq) and hgvs utils.
        exac :: annotation acquired via ExAC and dbsnp.

    Each component represents the most deleterious of the locus alleles.
    """
    var: SimpleVariant
    vcf: VcfVariantAnnotation
    hgvs: HgvsVariantAnnotation
    exac: ExacVariantAnnotation

    @classmethod
    def from_vcf_locus(cls, hgvs_machinery: HgvsMachinery, locus: _Record) -> 'VariantAnnotation':
        """
        :param hgvs_machinery: instance of hgvs machinery.
        :param locus: vcf record for locus.
        :return: variant annotation.
        """
        # generate simple variants for each allele
        variants = simple_variants_from_record(locus)

        # use hgvs to assess impact of each allele | pick the most deleterious allele
        variant, hgvs_ann = max(
            ((variant, HgvsVariantAnnotation.from_simple_variant(hgvs_machinery, variant)) for variant in variants),
            key=lambda var_hgvs: var_hgvs[1].feature_variant.impact if var_hgvs[1].feature_variant else -1)

        # use a 5'-normalized variant for ExAC
        variant_exac = hgvs_machinery.simple_variant_from_hgvs(hgvs_machinery.normalizer_5p.normalize(hgvs_ann.hgvs_g))

        return cls(
            var=variant,
            vcf=VcfVariantAnnotation.from_vcf_locus(locus, variant.alt_index),
            hgvs=hgvs_ann,
            exac=ExacVariantAnnotation.from_simple_variant(variant_exac))


def annotate_vcf(vcf_path: Path) -> Iterable[VariantAnnotation]:
    """
    Annotates a given VCF file.

    :param vcf_path: path of the vcf file.
    :return: variant annotations.
    """
    with vcf_path.open() as vcf_io:
        reader = Reader(vcf_io)
        hgvs_machinery = HgvsMachinery.from_assembly(assembly=TEMPUS_REFERENCE__ASSEMBLY[reader.metadata['reference']])
        yield from (VariantAnnotation.from_vcf_locus(hgvs_machinery, record) for record in reader)


def write_annotations_to_csv(annotations: Iterable[VariantAnnotation], out_io: TextIO):
    """
    Writes variant annotations to CSV.

    :param annotations: variant annotations
    :param out_io: open stream.
    """

    def coerce(value: Any) -> Any:
        if isinstance(value, (SequenceAlteration, FeatureVariant, SequenceVariant)):
            return str(value)
        elif isinstance(value, float):
            return round(value, 2)
        else:
            return value

    def to_dict(ann: VariantAnnotation) -> Dict[str, Any]:
        flat = flatten(asdict(ann), reducer=lambda key_a, key_b: key_b if key_a is None else '_'.join((key_a, key_b)))
        return {key: coerce(value) for key, value in flat.items()}

    writer = csv.DictWriter(out_io, fieldnames=(
        # variant fields
        'var_contig', 'var_pos', 'var_ref', 'var_alt', 'var_alt_index',
        # vcf fields
        'vcf_read_depth_site', 'vcf_read_depth_alt', 'vcf_perc_reads_alt', 'vcf_containing_samples',
        # hgvs fields
        'hgvs_sequence_alteration', 'hgvs_feature_variant', 'hgvs_hgvs_g', 'hgvs_hgvs_c', 'hgvs_hgvs_p', 'hgvs_gene',
        # exac fields
        'exac_allele_frequency', 'exac_consequences'))
    writer.writeheader()
    writer.writerows(map(to_dict, annotations))
