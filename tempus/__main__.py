"""
Usage:
    tempus annotate <VCF> <CSV> [--workers=<N>] [--chunk-size=<N>]
    tempus -h | --help

Arguments:
    VCF                 input variant call format file
    CSV                 output comma-separated file

Options:
    -h --help           show help
    --workers=<N>       number of workers [default: 1]
    --chunk-size=<N>    number of vcf records to process by each worker [default: 1]
"""
from concurrent.futures.process import ProcessPoolExecutor
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable

from docopt import docopt

from tempus.annotation import write_annotations_to_csv, annotate_vcf
from tempus.vcf import split_vcf


def _annotate_vcf(in_vcf: Path):
    annotations = annotate_vcf(in_vcf)
    with NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as out_csv_io:
        write_annotations_to_csv(annotations, out_csv_io)
    return Path(out_csv_io.name)


def _merge_csv(in_csvs: Iterable[Path], out_csv: Path):
    with out_csv.open('w') as out_io:
        for idx, in_csv in enumerate(in_csvs):
            with in_csv.open() as in_io:
                next(in_io) if idx > 0 else None
                out_io.writelines(in_io)
            in_csv.unlink()


def _handle_annotate(in_vcf: Path, out_csv: Path, max_workers: int, chunk_size: int):
    vcf_splits = split_vcf(in_vcf, chunk_size)
    with ProcessPoolExecutor(max_workers=max_workers) as proc_pool:
        tmp_csvs = proc_pool.map(_annotate_vcf, vcf_splits)
    _merge_csv(tmp_csvs, out_csv)


def cli():
    args = docopt(__doc__)
    if args['annotate']:
        _handle_annotate(
            in_vcf=Path(args['<VCF>']),
            out_csv=Path(args['<CSV>']),
            max_workers=int(args['--workers']),
            chunk_size=int(args['--chunk-size']))


if __name__ == '__main__':
    cli()
