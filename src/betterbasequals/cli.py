"""Module that contains the command line application."""

import argparse
from betterbasequals.get_good_bad_kmers import get_good_and_bad_kmers
from betterbasequals.utils import matches
from betterbasequals import __version__
import collections
from kmerpapa.algorithms import greedy_penalty_plus_pseudo
import numpy as np
from scipy.optimize import minimize


def greedy_papa_wrapper():
    def f(x):
        this_alpha, this_penalty, test_score = greedy_penalty_plus_pseudo.greedy_partition_CV(super_pattern, contextD, pseudo_counts, pargs, n_good, n_bad, penalty_values)
        return test_score[0]

def greedy_papa_NealderMead():
    x0 = np.array([1.0, 0.7, 0.8, 1.9, 1.2])
    res = minimize(rosen, x0, method='nelder-mead',
               options={'xatol': 1e-8, 'disp': True})

def get_parser():
    """
    Return the CLI argument parser.

    Returns:
        An argparse parser.
    """
    parser = argparse.ArgumentParser(
        prog="BetterBaseQuals",
        description='''
        Calculates sample-specific base qualities using overlapping reads.
        ''')

    parser.add_argument('-V', '--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument("bam_file", help="bam file")
    parser.add_argument("filter_bam_file", help="bam file from blood")
    parser.add_argument("twobit_file", help="Reference genome in two-bit format")
    parser.add_argument("--output_file_good", type=argparse.FileType('w'))
    parser.add_argument("--output_file_bad", type=argparse.FileType('w'))
    parser.add_argument("--radius", type=int, default=3)
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
        opts.radius,
    )

    super_pattern = 'N'*opts.radius + 'M' + 'N'*opts.radius
    if not opts.output_file_good is None:
        for kmer in  matches(super_pattern):
            print(kmer, good_kmers[kmer] , file = opts.output_file_good)
        opts.output_file_good.close()
    if not opts.output_file_bad is None:
        for kmer in matches(super_pattern):
            print(kmer, bad_kmers[kmer], file = opts.output_file_bad)
        opts.output_file_bad.close()

    pseudo_counts = [1,2,5,10,20,30,50,100,500,1000]
    penalty_values = range(1,15)
    n_good = sum(good_kmers.values())
    n_bad = sum(bad_kmers.values())
    print(n_good, n_bad)
    contextD = dict((x, (good_kmers[x], bad_kmers[x])) for x in matches(super_pattern))
    PapaArgs = collections.namedtuple('PapaArg', ('nfolds', 'iterations', 'seed',  'verbosity', 'CVfile'))
    pargs = PapaArgs(nfolds=10, iterations=5, seed=None, verbosity=1, CVfile=None)
    this_alpha, this_penalty, test_score = greedy_penalty_plus_pseudo.greedy_partition_CV(super_pattern, contextD, pseudo_counts, pargs, n_good, n_bad, penalty_values)
    print(this_alpha, this_penalty, test_score)
    return 0


