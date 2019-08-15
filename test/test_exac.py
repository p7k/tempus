import pytest

from tempus import SimpleVariant
from tempus.exac import annotate_variant_exac


@pytest.mark.parametrize('variant', [
    SimpleVariant(contig='1', pos=935222, ref='C', alt='A'),
], ids=['snp'])
def test_annotate_variant_exac(variant: SimpleVariant):
    ann = annotate_variant_exac(variant)
    assert ann.allele_frequency
    assert ann.consequences
