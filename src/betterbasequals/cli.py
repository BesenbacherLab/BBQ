"""Module that contains the command line application."""

import argparse
from betterbasequals.get_good_bad_kmers import get_good_and_bad_kmers, get_good_and_bad_kmers_w_filter
#from betterbasequals.call_mutations import MutationCaller
from betterbasequals.call_mutations import MutationValidator, BaseAdjuster
from betterbasequals.utils import *
from betterbasequals import __version__
from kmerpapa.algorithms import greedy_penalty_plus_pseudo
from kmerpapa.pattern_utils import get_M_U
from math import log10

# TODO: Har sat som constant nu. Lav flexible løsning.
VALID_BASEQUALS = [11,25,37]

def get_parser():
    """
    Return the CLI argument parser.

    Returns:
        An argparse parser.
    """
    parser = argparse.ArgumentParser(
        prog="bbq",
        description='''
        Calculates sample-specific base qualities using overlapping reads.
        ''')

    # top level args:
    parser.add_argument("--verbosity", type=int, default=1)
    parser.add_argument('-V', '--version', action='version', version=f'%(prog)s {__version__}')

    subparsers = parser.add_subparsers(dest='command', #help='commands',
        title="required commands",
        help='Select one of:')

    # args bam file:
    bam_parent = argparse.ArgumentParser(add_help=False)
    bam_parent.add_argument("--bam_file", required=True,
        help="bam file")
    bam_parent.add_argument("--twobit_file", required=True,
        help="Reference genome in two-bit format")
    bam_parent.add_argument('--region', '-r', type=str,
        help='only consider variants in this region')

    # args for filter bam file
    filter_parent = argparse.ArgumentParser(add_help=False)
    filter_parent.add_argument("--filter_bam_file", help="bam file from blood")

    # args for counting kmers:
    count_parent = argparse.ArgumentParser(add_help=False)
    count_parent.add_argument("--output_file_good", type=argparse.FileType('w'))
    count_parent.add_argument("--output_file_bad", type=argparse.FileType('w'))
    count_parent.add_argument("--radius", type=int, default=3)
    count_parent.add_argument('--min_depth', type=int, default=1,
        help="mminimum depth at a site to be considered as training data")
    count_parent.add_argument('--max_depth', type=int, default=5000,
        help="maximum depth at a site to be considered as training data")
    count_parent.add_argument('--filter_min_depth', type=int, default=1,
        help="minimum depth in filter_bam_file for a site to be considered as training data")
    count_parent.add_argument('--filter_max_depth', type=int, default=5000,
        help="maximum depth om filter_bam_file at a site to be considered as training data")


    # args for training models:    
    train_parent = argparse.ArgumentParser(add_help=False)
    train_parent.add_argument('--kmerpapa_method', type=str, default = "greedy",
        help='only consider variants in this region')
    train_parent.add_argument("--output_file_kmerpapa", type=argparse.FileType('w'))
    train_parent.add_argument('-N', '--nfolds', type=int, metavar='N', default=2,
        help='Number of folds to use when fitting hyperparameters in kmerpapa')
    train_parent.add_argument('-i', '--iterations', type=int, default=1, metavar='i',
        help='Repeat cross validation i times when fitting hyperparameters in kmerpapa')
    train_parent.add_argument('--seed', type=int,
        help='seed for numpy.random')

    # args for validating models:    
    validate_parent = argparse.ArgumentParser(add_help=False)
    validate_parent.add_argument("--validation_bam_file", help="hifi bam file")

    # args for printing polished bam:    
    adjust_parent = argparse.ArgumentParser(add_help=False)
    adjust_parent.add_argument('--outbam', type=str,
        help="Bam file with adjusted base qualities.", required=True)

    # args for calling somatic variants:    
    call_parent = argparse.ArgumentParser(add_help=False)
    call_parent.add_argument('--outfile', type=str,
        help="output file")

    count_parser = subparsers.add_parser('count', 
        description='Count good and bad k-mers',
        help =  'Count good and bad k-mers',
        parents=[bam_parent, filter_parent, count_parent],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    train_parser = subparsers.add_parser('train', 
        description='First run "count" then train model to distinguish good and bad k-mers.', 
        help = 'First run "count" then train model to distinguish good and bad k-mers.',
        parents=[bam_parent, filter_parent, count_parent, train_parent])

    validate_parser = subparsers.add_parser('validate', 
        description = 'First run "count" and "train" then print validation data',
        help = 'First run "count" and "train" then print validation data',
        parents=[bam_parent, filter_parent, count_parent, train_parent, validate_parent])

    adjust_parser = subparsers.add_parser('adjust', 
        description = 'First run "count" and "train" then output bam with adjusted base qualities.', 
        help = 'First run "count" and "train" then output bam with adjusted base qualities.', 
        parents = [bam_parent, filter_parent, count_parent, train_parent, adjust_parent])

    call_parser = subparsers.add_parser('call', 
        description = 'First run "count" and "train" then call variants', 
        help = 'First run "count" and "train" then call variants', 
        parents = [bam_parent, count_parent, train_parent, call_parent])

    train_only_parser = subparsers.add_parser('train_only', 
        description = 'Train model to distinguish good and bad k-mers.',
        help = 'Train model to distinguish good and bad k-mers.',
        parents = [train_parent])
    train_only_parser.add_argument("--input_file_good", type=argparse.FileType('r'))
    train_only_parser.add_argument("--input_file_bad", type=argparse.FileType('r'))

    validate_only_parser = subparsers.add_parser('validate_only', 
        description = 'Print validation data.',
        help = 'Print validation data.', 
        parents = [bam_parent, filter_parent, validate_parent]) 
    validate_only_parser.add_argument("--input_file_kmerpapa", type=argparse.FileType('r'))
    
    adjust_only_parser = subparsers.add_parser('adjust_only', 
        description = 'Output bam with adjusted base qualities.', 
        help = 'Output bam with adjusted base qualities.',
        parents = [bam_parent, adjust_parent])
    adjust_only_parser.add_argument("--input_file_kmerpapa", type=argparse.FileType('r'))

    call_only_parser = subparsers.add_parser('call_only', 
        description = 'Call variants',
        help = 'Call variants',
        parents=[call_parent])

    return parser


def run_get_good_and_bad(opts):
    if opts.verbosity > 0:
        eprint("Counting good and bad kmers")
    good_kmers, bad_kmers = get_good_and_bad_kmers(
        opts.bam_file,
        opts.twobit_file,
        opts.chrom,
        opts.start,
        opts.end,
        opts.min_depth,
        opts.max_depth,
        opts.radius,
    )
    print_good_and_bad(opts, good_kmers, bad_kmers)
    return good_kmers, bad_kmers

def run_get_good_and_bad_w_filter(opts):
    if opts.verbosity > 0:
        eprint("Counting good and bad kmers")
    good_kmers, bad_kmers = get_good_and_bad_kmers_w_filter(
        opts.bam_file,
        opts.twobit_file,
        opts.filter_bam_file,
        VALID_BASEQUALS,
        opts.chrom,
        opts.start,
        opts.end,
        opts.min_depth,
        opts.max_depth,
        opts.radius,
    )
    print_good_and_bad(opts, good_kmers, bad_kmers)
    return good_kmers, bad_kmers


def run_get_kmerpapas(opts, good_kmers, bad_kmers):
    if opts.verbosity > 0:
        eprint("Training kmer pattern partitions")
    kmer_papas = {}
    #pseudo_counts = [1,2,5,10,20,30,50,100,500,1000]
    #penalty_values = range(1,15)
    for bqual in good_kmers:
        kmer_papas[bqual] = {}
        eprint(f'Handling base_qual: {bqual}')
        for mtype in good_kmers[bqual]:
            eprint(f'Handling mutation type: {mtype}')
            radius = len(next(iter(good_kmers[bqual]["A->C"].keys())))//2
            super_pattern = 'N'*radius + mtype[0] + 'N'*radius
            n_good = sum(good_kmers[bqual][mtype].values())
            n_bad = sum(bad_kmers[bqual][mtype].values())
            eprint(mtype, n_bad, n_good)
            contextD = dict((x, (bad_kmers[bqual][mtype][x], good_kmers[bqual][mtype][x])) for x in matches(super_pattern))
            CV = greedy_penalty_plus_pseudo.BaysianOptimizationCV(super_pattern, contextD, opts.nfolds, opts.iterations, opts.seed)
            best_alpha, best_penalty, test_score = CV.get_best_a_c()
            my=n_bad/(n_good+n_bad)
            best_beta = (best_alpha*(1.0-my))/my
            eprint(best_penalty, best_alpha, best_beta, my)
            best_score, M, U, names = greedy_penalty_plus_pseudo.greedy_partition(super_pattern, contextD, best_alpha, best_beta, best_penalty, {})
            counts = []
            for pat in names:
                counts.append(get_M_U(pat, contextD))
            kmer_papas[bqual][mtype] = {}
            for i in range(len(names)):
                pat = names[i]
                M, U = counts[i]
                p = (M + best_alpha)/(M + U + best_alpha + best_beta)
                if opts.verbosity > 0:
                    eprint(mtype, pat, p, -10*log10(p/(1-p)))
                if not opts.output_file_kmerpapa is None:
                    print(bqual, mtype, pat, -10*log10(p), file=opts.output_file_kmerpapa)
                    #print(bqual, mtype, pat, -10*log10(p/(1-p)), file=opts.output_file_kmerpapa)
                for context in matches(pat):
                    kmer_papas[bqual][mtype][context] = -10*log10(p)
    if not opts.output_file_kmerpapa is None:
        opts.output_file_kmerpapa.close()
    return kmer_papas

def run_validation(opts, kmer_papas):
    if opts.verbosity > 0:
        eprint("Printing validation data")
    validator = \
        MutationValidator(
            opts.bam_file, 
            opts.filter_bam_file, 
            opts.validation_bam_file, 
            opts.twobit_file, 
            kmer_papas)
    
    validator.call_mutations(opts.chrom, opts.start, opts.end)

def run_adjust(opts, kmer_papas):
    if opts.verbosity > 0:
        eprint("Adjusting base qualities")
    adjuster = \
        BaseAdjuster(
            opts.bam_file,
            opts.twobit_file, 
            kmer_papas,
            opts.outbam)
    
    n_corrections, n_corrected_reads, n_uncorrected_reads, n_filtered = \
        adjuster.call_mutations(opts.chrom, opts.start, opts.end)

    if opts.verbosity > 0:
        eprint(f'corrected {n_corrections} base qualities in {n_corrected_reads} reads')
        eprint(f'{n_uncorrected_reads} reads were written with no corrections')
        eprint(f'{n_filtered} reads were filtered')


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

    if "region" in opts:
        parse_opts_region(opts)

    if opts.command is None:
        parser.print_help(sys.stderr)
        return 1


    if not opts.command in ['train_only', 'validate_only', 'call_only', 'adjust_only']:
        #if opts.filter_bam_file is None:            
        #    good_kmers, bad_kmers = run_get_good_and_bad(opts)
        #else:
        good_kmers, bad_kmers = run_get_good_and_bad_w_filter(opts)
    elif opts.command == 'train_only':
        good_kmers, bad_kmers = read_kmers(opts)

    if opts.command == 'count':
        return 0

    if not opts.command in ['validate_only', 'call_only', 'adjust_only']:
        kmer_papas = run_get_kmerpapas(opts, good_kmers, bad_kmers)
    else:
        eprint("Reading kmer pattern partitions")
        kmer_papas = read_kmer_papas(opts)

    # for x in kmer_papas:
    #     eprint(x)

    if opts.command in ['train_only', 'train']:
        return 0

    if opts.command in ['validate', 'validate_only']:
        run_validation(opts, kmer_papas)
    elif opts.command in ['call', 'call_only']:
        eprint("Call not implemented yet")
    elif opts.command in ['adjust', 'adjust_only']:
        run_adjust(opts, kmer_papas)
    else:
        eprint("Unknown command: {opts.command}")
        return 1

    #caller = MutationCaller(opts.bam_file, opts.filter_bam_file, opts.twobit_file, kmer_papas)

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


