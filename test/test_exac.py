import pytest

from tempus import SimpleVariant
from tempus.exac import ExacVariantAnnotation


@pytest.mark.parametrize('variant', [
    SimpleVariant(contig='1', pos=935222, ref='C', alt='A'),
], ids=['snp'])
def test_annotate_variant_exac(variant: SimpleVariant):
    ann = ExacVariantAnnotation.from_simple_variant(variant)
    assert ann.allele_frequency
    assert ann.consequences
