"""Module that contains the command line application."""

import argparse
from get_good_bad_kmers import get_good_and_bad_kmers

def get_parser():
    """
    Return the CLI argument parser.

    Returns:
        An argparse parser.
    """
    parser = argparse.ArgumentParser(
        prog="BetterBaseQuals",
        description='''
        
        ''')
    parser.add_argument('-V', '--version', action='store_true',
        help='Print version number and return')
    parser.add_argument("bam_file", help="bam file")
    parser.add_argument("filter_bam_file", help="bam file from blood")
    parser.add_argument("twobit_file", help="Reference genome in two-bit format")
    parser.add_argument("output_file_good")
    parser.add_argument("output_file_bad")
    parser.add_argument('--region','-r',type=str,help='only consider variants in this region')
    return parser


def main(args = None):
    """
    Run the main program.

    This function is executed when you type `BetterBaseQuals` or `python -m betterbasequals`.

    Arguments:
        args: Arguments passed from the command line.

    Returns:
        An exit code.
    """
    parser = get_parser()
    opts = parser.parse_args(args=args)
    print(opts) 

    if opts.version:
        from betterbasequals import __version__
        print("version:", __version__)
        print()
        return 0

    if opts.region is None:
        chrom = None
        start = None
        end = None
    else:
        chrom, end_points = opts.region.split(':')
        start, end = end_points.split('-')
        start = int(start)
        end = int(end)

    good_kmers, bad_kmers = get_good_and_bad_kmers(
        opts.bam_file,
        opts.filter_bam_file,
        opts.twobit_file,
        chrom,
        start,
        end,
    )

    with open(opts.output_file_good, 'w') as f:
        for kmer in matches("NNNMNNN"):
            print(kmer, good_kmers[kmer], file = f)
    with open(opts.output_file_bad, 'w') as f:
        for kmer in matches("NNNMNNN"):
            print(kmer, bad_kmers[kmer], file = f)
    
    return 0
