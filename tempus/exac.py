"""
ExAC-based annotation functions - using the non-bulk REST API.
"""
from dataclasses import dataclass
from typing import Optional, List

import requests

from tempus import SimpleVariant

API_BASE_URL = 'http://exac.hms.harvard.edu/rest/'


@dataclass(frozen=True)
class ExacVariantAnnotation:
    allele_frequency: Optional[float]
    consequences: List[str]

    @classmethod
    def from_simple_variant(cls, variant: SimpleVariant) -> 'ExacVariantAnnotation':
        """
        Annotate variant using ExAC REST API.

        :param variant: 5p-normalized (left shuffled) simple variant. [base coordinate system]
        :return: variant annotation
        """
        var = f'{variant.contig}-{variant.pos}-{variant.ref}-{variant.alt}'
        response = requests.get(url=f'{API_BASE_URL}/variant/{var}')
        response.raise_for_status()
        data = response.json()

        return cls(
            allele_frequency=data.get('variant', {}).get('allele_freq'),
            consequences=list((data.get('consequence') or {}).keys()))
