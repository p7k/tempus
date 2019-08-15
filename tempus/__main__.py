"""
Usage:
    tempus annotate VCF [-o <CSV>]

Arguments:
    VCF         vcf file

Options:
    -h --help   show help
    -o CSV      output csv
"""
import sys
from pathlib import Path
from typing import Any, Dict

from docopt import docopt

from tempus.annotation import write_annotations_to_csv, annotate_vcf


def handle_annotate(args: Dict[str, Any]):
    annotations = annotate_vcf(Path(args['VCF']))
    output = args.get('o', '-')
    if output != '-':
        with open(output, 'w', newline='') as io_out:
            write_annotations_to_csv(annotations, io_out)
    else:
        write_annotations_to_csv(annotations, sys.stdout)


def cli():
    args = docopt(__doc__)

    try:
        if args['annotate']:
            handle_annotate(args)
    except KeyboardInterrupt:
        pass


if __name__ == '__main__':
    cli()
