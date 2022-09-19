"""Module that contains the command line application."""

import argparse
from betterbasequals.get_good_bad_kmers import get_good_and_bad_kmers
#from betterbasequals.call_mutations import MutationCaller
from betterbasequals.call_mutations import MutationValidator
from betterbasequals.utils import matches, mtypes, eprint
from betterbasequals import __version__
from kmerpapa.algorithms import greedy_penalty_plus_pseudo
from kmerpapa.pattern_utils import get_M_U
from math import log10


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
    parser.add_argument("validation_bam_file", help="hifi bam file")
    parser.add_argument("twobit_file", help="Reference genome in two-bit format")
    parser.add_argument("--output_file_good", type=argparse.FileType('w'))
    parser.add_argument("--output_file_bad", type=argparse.FileType('w'))
    parser.add_argument("--output_file_kmerpapa", type=argparse.FileType('w'))
    parser.add_argument("--radius", type=int, default=3)
    parser.add_argument('--min_depth', type=int, default=1,
        help="mminimum depth at a site to be considered as training data")
    parser.add_argument('--max_depth', type=int, default=1,
        help="maximum depth at a site to be considered as training data")
    parser.add_argument('--region', '-r', type=str,
        help='only consider variants in this region')
    parser.add_argument('--outbam', type=str,
        help="Bam file with adjusted base qualities.")
    parser.add_argument('-N', '--nfolds', type=int, metavar='N', default=2,
        help='Number of folds to use when fitting hyperparameters in kmerpapa')
    parser.add_argument('-i', '--iterations', type=int, default=1, metavar='i',
        help='Repeat cross validation i times when fitting hyperparameters in kmerpapa')
    parser.add_argument('--seed', type=int,
        help='seed for numpy.random')
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

    if opts.region is None:
        chrom = None
        start = None
        end = None
    else:
        chrom, end_points = opts.region.split(':')
        start, end = end_points.split('-')
        start = int(start)
        end = int(end)

    #outfile = pysam.AlignmentFile(opts.outbam, "w", template=opts.bam_file)
    eprint("Count good and bad kmers")
    good_kmers, bad_kmers = get_good_and_bad_kmers(
        opts.bam_file,
        opts.twobit_file,
        chrom,
        start,
        end,
        opts.min_depth,
        opts.max_depth,
        opts.radius,
    )

    if not opts.output_file_good is None:
        for mtype in mtypes:
            super_pattern = 'N'*opts.radius + mtype[0] + 'N'*opts.radius
            for kmer in  matches(super_pattern):
                print(mtype, kmer, good_kmers[mtype][kmer] , file = opts.output_file_good)
        opts.output_file_good.close()
    if not opts.output_file_bad is None:
        for mtype in mtypes:
            super_pattern = 'N'*opts.radius + mtype[0] + 'N'*opts.radius
            for kmer in matches(super_pattern):
                print(mtype, kmer, bad_kmers[mtype][kmer], file = opts.output_file_bad)
        opts.output_file_bad.close()


    kmer_papas = {}
    #pseudo_counts = [1,2,5,10,20,30,50,100,500,1000]
    #penalty_values = range(1,15)
    for mtype in mtypes:
        eprint(f'Handling {mtype}')
        super_pattern = 'N'*opts.radius + mtype[0] + 'N'*opts.radius
        n_good = sum(good_kmers[mtype].values())
        n_bad = sum(bad_kmers[mtype].values())
        print(mtype, n_bad, n_good)
        contextD = dict((x, (bad_kmers[mtype][x], good_kmers[mtype][x])) for x in matches(super_pattern))
        CV = greedy_penalty_plus_pseudo.BaysianOptimizationCV(super_pattern, contextD, opts.nfolds, opts.iterations, opts.seed)
        best_alpha, best_penalty, test_score = CV.get_best_a_c()
        my=n_bad/(n_good+n_bad)
        best_beta = (best_alpha*(1.0-my))/my
        print(best_penalty, best_alpha, best_beta, my)
        best_score, M, U, names = greedy_penalty_plus_pseudo.greedy_partition(super_pattern, contextD, best_alpha, best_beta, best_penalty, {})
        counts = []
        for pat in names:
            counts.append(get_M_U(pat, contextD))
        kmer_papas[mtype] = {}
        for i in range(len(names)):
            pat = names[i]
            M, U = counts[i]
            p = (M + best_alpha)/(M + U + best_alpha + best_beta)
            if not opts.output_file_kmerpapa is None:
                print(mtype, pat, p, -10*log10(p/(1-p)))
                print(mtype, pat, p, -10*log10(p/(1-p)), file=opts.output_file_kmerpapa)
            for context in matches(pat):
                kmer_papas[mtype][context] = -10*log10(p/(1-p))
    #caller = MutationCaller(opts.bam_file, opts.filter_bam_file, opts.twobit_file, kmer_papas)

    validator = \
        MutationValidator(
            opts.bam_file, 
            opts.filter_bam_file, 
            opts.validation_bam_file, 
            opts.twobit_file, 
            kmer_papas)
    
    validator.call_mutations(chrom, start, end)


    # run_mutation_caller(
    #     opts.bam_file, 
    #     opts.filter_bam_file, 
    #     opts.twobit_file, 
    #     chrom, 
    #     start, 
    #     end, 
    #     opts.radius,
    # )

    return 0


