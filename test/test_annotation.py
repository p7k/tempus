from pathlib import Path

import pytest

from tempus.annotation import annotate_vcf


@pytest.mark.parametrize('vcf_path', [Path('test/data/data_sys_test.vcf')])
def test_an(vcf_path: Path):
    anns = tuple(annotate_vcf(vcf_path))
    assert anns
