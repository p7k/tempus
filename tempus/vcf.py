"""
Convenience functions for working with the `jamescasbon/PyVCF` library.
"""
from dataclasses import dataclass
from typing import Iterable, Tuple, List, Dict

from vcf.model import _Record, _Call

from tempus import SimpleVariant, Assembly

TEMPUS_REFERENCE__ASSEMBLY: Dict[str, Assembly] = {
    '/data/human_g1k_v37.fasta': Assembly.GRCH37,
}


def simple_variants_from_record(record: _Record) -> Iterable[Tuple[int, SimpleVariant]]:
    """
    Produces SimpleVariants for each ALT allele at a locus.

    :param record: vcf
    :return: allele index and SimpleVariant
    """
    return (
        (idx, SimpleVariant(contig=record.CHROM, pos=record.POS, ref=record.REF, alt=alt.sequence))
        for idx, alt in enumerate(record.ALT, 1))


def calls_containing_allele(record: _Record, allele_index: int) -> Iterable[_Call]:
    """
    :param record: vcf site object.
    :param allele_index: zero-based allele index.
    :return:
    """
    return (
        call for call in record.samples
        if any(int(gt_allele) == allele_index for gt_allele in call.gt_alleles))


def read_depth_allele(record: _Record, allele_index: int) -> int:
    """
    :param record: site object.
    :param allele_index: zero-based allele index
    :return: sum of allele read depths across all samples.
    """
    return sum(call.data.DPR[allele_index] for call in record.samples)


@dataclass(frozen=True)
class VcfVariantAnnotation:
    chrom: str
    pos: int
    ref: str
    alt: List[str]
    allele_index: int
    read_depth_site: int
    read_depth_alt: int
    perc_reads_alt: float
    containing_samples: List[str]


def annotate_variant_vcf(record: _Record, allele_index: int) -> VcfVariantAnnotation:
    read_depth_ref = read_depth_allele(record, 0)
    read_depth_alt = read_depth_allele(record, allele_index)

    return VcfVariantAnnotation(
        chrom=record.CHROM,
        pos=record.POS,
        ref=record.REF,
        alt=list(map(str, record.ALT)),
        allele_index=allele_index,
        read_depth_site=record.INFO['DP'],
        read_depth_alt=read_depth_alt,
        perc_reads_alt=float(read_depth_alt) / read_depth_ref,
        containing_samples=[call.sample for call in calls_containing_allele(record, allele_index)])
