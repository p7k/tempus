import csv
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, Tuple, Dict, Any, TextIO

from flatten_dict import flatten
from hgvs.sequencevariant import SequenceVariant
from vcf import Reader
from vcf.model import _Record

from tempus import SimpleVariant, FeatureVariant, SequenceAlteration, Assembly
from tempus.exac import ExacVariantAnnotation, annotate_variant_exac
from tempus.hgvs import annotate_variant_hgvs, HgvsVariantAnnotation, HgvsMachinery
from tempus.vcf import simple_variants_from_record, VcfVariantAnnotation, annotate_variant_vcf, assembly_from_vcf


@dataclass(frozen=True)
class VariantAnnotation:
    vcf: VcfVariantAnnotation
    hgvs: HgvsVariantAnnotation
    exac: ExacVariantAnnotation


def annotate_vcf_record(hgvs_machinery: HgvsMachinery, record: _Record) -> VariantAnnotation:
    # use hgvs to assess impact
    hgvs_anns: Iterable[Tuple[int, SimpleVariant, HgvsVariantAnnotation]] = (
        (allele_idx, variant, annotate_variant_hgvs(hgvs_machinery, variant))
        for allele_idx, variant in simple_variants_from_record(record))

    # pick most deleterious
    allele_idx, variant, hgvs_ann = max(
        hgvs_anns, key=lambda idx_ann: idx_ann[2].feature_variant.impact if idx_ann[2].feature_variant else -1)

    # prepare a 5p-normalized variant for ExAC
    variant_exac = hgvs_machinery.simple_variant_from_hgvs(hgvs_machinery.normalizer_5p.normalize(hgvs_ann.hgvs_g))

    return VariantAnnotation(
        vcf=annotate_variant_vcf(record, allele_idx),
        hgvs=hgvs_ann,
        exac=annotate_variant_exac(variant_exac))


def _proc_init(assembly: Assembly):
    globals()['HGVS_MACHINERY'] = HgvsMachinery.from_assembly(assembly)


def _task_wrapper(record: _Record) -> VariantAnnotation:
    hgvs_machinery = globals()['HGVS_MACHINERY']
    return annotate_vcf_record(hgvs_machinery, record)


import pickle


def annotate_vcf(vcf_path: Path) -> Iterable[VariantAnnotation]:
    assembly = assembly_from_vcf(vcf_path)

    with vcf_path.open() as vcf_io:
        records = pickle.loads(pickle.dumps(tuple(Reader(vcf_io))))

    _proc_init(assembly)
    yield from map(_task_wrapper, records)

    #
    # with ProcessPoolExecutor(max_workers=1, initializer=_proc_init, initargs=[assembly]) as pool_exec:
    #     yield from pool_exec.map(_task_wrapper, records, chunksize=1000, timeout=3)


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
        'vcf_chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt',
        'vcf_read_depth_site', 'vcf_read_depth_alt', 'vcf_perc_reads_alt', 'vcf_containing_samples',
        # hgvs fields
        'hgvs_sequence_alteration', 'hgvs_feature_variant', 'hgvs_hgvs_g', 'hgvs_hgvs_c', 'hgvs_hgvs_p', 'hgvs_gene',
        # exac fields
        'exac_allele_frequency', 'exac_consequences'))
    writer.writeheader()
    writer.writerows(map(to_dict, annotations))
