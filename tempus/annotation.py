import csv
from concurrent.futures.process import ProcessPoolExecutor
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, Tuple, Dict, Any, TextIO, Sequence

from flatten_dict import flatten
from hgvs.sequencevariant import SequenceVariant
from vcf import Reader

from tempus import SimpleVariant, FeatureVariant, SequenceAlteration, Assembly
from tempus.exac import ExacVariantAnnotation
from tempus.hgvs import HgvsVariantAnnotation, HgvsMachinery
from tempus.vcf import simple_variants_from_record, TEMPUS_REFERENCE__ASSEMBLY, VcfVariantAnnotation, \
    annotate_variant_vcf


@dataclass(frozen=True)
class VariantAnnotation:
    """
    Locus-level annotation encapsulating the following:

        vcf :: annotations extracted directly from the vcf itself
        hgvs :: annotations acquired via UTA (RefSeq) and hgvs utils.
        exac :: annotation acquired via ExAC and dbsnp.

    Each component represents the most deleterious of the site's alleles.
    """
    vcf: VcfVariantAnnotation
    hgvs: HgvsVariantAnnotation
    exac: ExacVariantAnnotation

    # TODO preparing VcfAnnotation(s) for every allele is a workaround for pickl(ing) limitations of pyvcf.Record
    # TODO it's ugly and should be refactored
    @classmethod
    def from_simple_variants(
            cls, hgvs_machinery: HgvsMachinery,
            allele_variants: Iterable[Tuple[SimpleVariant, VcfVariantAnnotation]]) -> 'VariantAnnotation':
        """
        :param hgvs_machinery: instance of hgvs machinery.
        :param allele_variants: simple variants coupled with vcf annotations for each allele at a locus.
        :return: variant annotation.
        """
        # use hgvs to assess impact of each allele
        hgvs_anns = (
            (variant, vcf_ann, HgvsVariantAnnotation.from_simple_variant(hgvs_machinery, variant))
            for variant, vcf_ann in allele_variants)

        # pick most deleterious allele
        variant, vcf_ann, hgvs_ann = max(
            hgvs_anns, key=lambda idx_ann: idx_ann[2].feature_variant.impact if idx_ann[2].feature_variant else -1)

        # use a 5'-normalized variant for ExAC
        variant_exac = hgvs_machinery.simple_variant_from_hgvs(hgvs_machinery.normalizer_5p.normalize(hgvs_ann.hgvs_g))
        exac_ann = ExacVariantAnnotation.from_simple_variant(variant_exac)

        return cls(vcf_ann, hgvs_ann, exac_ann)


def _proc_init(assembly: Assembly):
    globals()['HGVS_MACHINERY'] = HgvsMachinery.from_assembly(assembly)


def _task_wrapper(allele_variants: Sequence[Tuple[SimpleVariant, VcfVariantAnnotation]]) -> VariantAnnotation:
    hgvs_machinery = globals()['HGVS_MACHINERY']
    return VariantAnnotation.from_simple_variants(hgvs_machinery, allele_variants)


def annotate_vcf(vcf_path: Path, max_workers=1) -> Iterable[VariantAnnotation]:
    """
    Annotates a given VCF file.

    :param vcf_path: path of the vcf file.
    :param max_workers: number of parallel processes.
    :return: variant annotations.
    """
    with vcf_path.open() as vcf_io:
        reader = Reader(vcf_io)
        assembly = TEMPUS_REFERENCE__ASSEMBLY[reader.metadata['reference']]

        # TODO thread pool would seem more appropriate if we are network IO bounded (requests to UTA, NCBI, and ExAC)
        # TODO but HGVS can work with a local seq repo `biocommons/biocommons.seqrepo`, but it's not thread-safe
        # TODO hence the proc pool, which might be a bit of a premature optimization and leads to messier APIs.
        with ProcessPoolExecutor(max_workers=max_workers, initializer=_proc_init, initargs=[assembly]) as pool_exec:
            # TODO see note on VariantAnnotation regarding over-engineered API
            allele_variants = (tuple(
                (simple_variant, annotate_variant_vcf(record, allele_idx))
                for allele_idx, simple_variant in simple_variants_from_record(record)) for record in reader)
            yield from pool_exec.map(_task_wrapper, allele_variants, chunksize=100, timeout=60)


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
        # vcf fields
        'vcf_chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'vcf_allele_index',
        'vcf_read_depth_site', 'vcf_read_depth_alt', 'vcf_perc_reads_alt', 'vcf_containing_samples',
        # hgvs fields
        'hgvs_sequence_alteration', 'hgvs_feature_variant', 'hgvs_hgvs_g', 'hgvs_hgvs_c', 'hgvs_hgvs_p', 'hgvs_gene',
        # exac fields
        'exac_allele_frequency', 'exac_consequences'))
    writer.writeheader()
    writer.writerows(map(to_dict, annotations))
