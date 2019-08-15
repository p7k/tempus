from pathlib import Path

import pytest
from vcf import Reader


@pytest.mark.parametrize('vcf_path', [Path('test/data/Challenge_data.vcf')])
def test_parse_vcf(vcf_path: Path):
    with vcf_path.open() as vcf_io:
        reader = Reader(vcf_io)
        for locus in reader:
            assert locus
