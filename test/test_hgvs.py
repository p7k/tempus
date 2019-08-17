import pytest

from tempus import Assembly, SimpleVariant
from tempus.hgvs import HgvsMachinery


@pytest.fixture(scope='module')
def hgvs_machinery() -> HgvsMachinery:
    return HgvsMachinery.from_assembly(Assembly.GRCH37)


@pytest.mark.parametrize('variant, hgvs_g_str', [
    pytest.param(SimpleVariant(contig='1', pos=100000, ref='C', alt='G'), 'NC_000001.10:g.100000C>G', id='snp'),
    pytest.param(SimpleVariant(contig='X', pos=100000, ref='C', alt='CG'), 'NC_000023.10:g.100000delinsCG', id='ins'),
    pytest.param(
        SimpleVariant(contig='Y', pos=100000, ref='CA', alt='C'), 'NC_000024.9:g.100000_100001delinsC', id='del'),
])
def test_hgvs_simple_variant_conversions(hgvs_machinery: HgvsMachinery, variant: SimpleVariant, hgvs_g_str: str):
    """
    Note: hgvs is not normalized.
    """
    hgvs_g = hgvs_machinery.hgvs_from_simple_variant(variant, vtype='g')
    assert str(hgvs_g) == hgvs_g_str

    _variant = hgvs_machinery.simple_variant_from_hgvs(hgvs_g)
    assert _variant == variant
