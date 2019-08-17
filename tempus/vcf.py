"""
Convenience functions for working with the `jamescasbon/PyVCF` library and annotation functions which draw from the VCF.
"""
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable, List, Dict

from more_itertools import partition, chunked
from vcf.model import _Record, _Call

from tempus import SimpleVariant, Assembly

TEMPUS_REFERENCE__ASSEMBLY: Dict[str, Assembly] = {
    '/data/human_g1k_v37.fasta': Assembly.GRCH37,
}


def simple_variants_from_record(record: _Record) -> Iterable[SimpleVariant]:
    """
    Produces SimpleVariants for each ALT allele at a locus.

    :param record: vcf
    :return: simple variants
    """
    return (
        SimpleVariant(contig=record.CHROM, pos=record.POS, ref=record.REF, alt=alt.sequence, alt_index=idx)
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
    read_depth_site: int
    read_depth_alt: int
    perc_reads_alt: float
    containing_samples: List[str]

    @classmethod
    def from_vcf_locus(cls, locus: _Record, allele_index: int) -> 'VcfVariantAnnotation':
        """
        Annotate a particular allele by extracting data directly from VCF.

        :param locus: vcf record for a given locus.
        :param allele_index: zero-based index for the allele to be annotated.
        :return: allele annotation.
        """
        read_depth_ref = read_depth_allele(locus, 0)
        read_depth_alt = read_depth_allele(locus, allele_index)
        return cls(
            read_depth_site=locus.INFO['DP'],
            read_depth_alt=read_depth_alt,
            perc_reads_alt=float(read_depth_alt) / read_depth_ref,
            containing_samples=[call.sample for call in calls_containing_allele(locus, allele_index)])


def split_vcf(vcf_path: Path, chunk_size: int) -> Iterable[Path]:
    """
    A simple utility for splitting a VCF file into chunk-sized VCF files.

    Note: all splits keep the original header.

    :param vcf_path: input vcf file.
    :param chunk_size: max number of records in each file.
    :return: paths of tmp split files.
    """
    with vcf_path.open() as vcf_io:
        records, header = partition(lambda line: line.startswith('#'), vcf_io)
        for n, chunk in enumerate(chunked(records, chunk_size)):
            header = header if isinstance(header, tuple) else tuple(header)
            with NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as out_io:
                out_io.writelines(header)
                out_io.writelines(chunk)
            yield Path(out_io.name)
