#! /usr/bin/env python

from argparse import ArgumentParser

from serutils.sequence.fastq_utils import subsample_SC



def main():
    opts = get_options()

    subsample_SC(fname=opts.fastq_file, outfname=opts.output, 
                 reading_chunk=opts.chunk_size, ncells=opts.ncells,
                 sub_sample_size=opts.sub_sample_size, regex=opts.regex)


def get_options():
    parser = ArgumentParser()

    parser.add_argument('-i', dest='fastq_file', metavar='PATH', required=True,
                        default=False, help='input single-cell FASTQ file')
    parser.add_argument('-o', '--output', dest='output', metavar='PATH', required=True,
                        help='''Output subsampled FASTQ file''')
    parser.add_argument('-s', '--sub_sample_size', dest='sub_sample_size', default=1_000, type=int,
                        help='''[%(default)s] number of entries wanted in the output FASTQ.''')
    parser.add_argument('--ncells', dest='ncells', default="96", type=int,
                        help='''[%(default)s] number of sequenced cells in the FASTQ.''')
    parser.add_argument('--chunk_size', dest='chunk_size', default=1_000, type=int,
                        help='''[%(default)s] number of byte to read at once 
                        (2 entries should be contained in this length).''')
    parser.add_argument('--regex', dest='regex', default="^(@).*:([0-9]+):[0-9]+:[0-9]+ .*\n[ACGTN]+\n\+\n.*$",
                        help='''[%(default)s] regexp pattern to extract cell-ID from FASTQ read-ID field''')
    opts = parser.parse_args()
    return opts

if __name__ == "__main__":
    exit(main())