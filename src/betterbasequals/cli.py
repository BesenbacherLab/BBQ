"""Module that contains the command line application."""

import argparse
from betterbasequals.get_good_bad_kmers import MutationCounterWFilter
from betterbasequals.model_validators import MutationValidator, ListMutationValidator
from betterbasequals.bam_adjusters import BaseAdjuster
from betterbasequals.somatic_callers import SomaticMutationCaller
from betterbasequals.utils import *
from betterbasequals import __version__
from betterbasequals.kmerpapa_utils import get_kmerpapa

# TODO: Possible input basequalities are hardcoded. Should make flexible solution.
#VALID_BASEQUALS = [11,25,37]

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
    
    read_filter_parent = argparse.ArgumentParser(add_help=False)
    read_filter_parent.add_argument('--min_enddist', type=int, default=6, metavar="M",
        help="Ignore bases in the first M or last M positions in the read")
    read_filter_parent.add_argument('--max_mismatch', type=int, default=1, metavar="M",
        help="Ignore alt reads if the read has more than M mismatches to the reference")

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
        help='algorithm to use for calculating kmer-papas',
        choices=['greedy', 'optimal'])
    train_parent.add_argument('--correction_type', type=str, default = "bad_vs_no",
        help='should we compare bad variants to "good variants"(SNVs) or to "no variant" (homozygous ref sites)',
        choices=["bad_vs_good", "bad_vs_no"])
    train_parent.add_argument("--output_file_kmerpapa", type=argparse.FileType('w'))
    train_parent.add_argument('-N', '--nfolds', type=int, metavar='N', default=2,
        help='Number of folds to use when fitting hyperparameters in kmerpapa')
    train_parent.add_argument('-i', '--iterations', type=int, default=1, metavar='i',
        help='Repeat cross validation i times when fitting hyperparameters in kmerpapa')
    train_parent.add_argument('--seed', type=int,
        help='seed for numpy.random')
    train_parent.add_argument('--same_good', action='store_true')
    train_parent.add_argument(
        '-a', '--pseudo_counts', type=float, metavar='a', nargs='+', default = [1,10],
        help='Different pseudo count (alpha) values to test using cross validation')
    train_parent.add_argument(
        '-c', '--penalty_values', type=float, metavar='c', nargs='+', default = [5,10,15,20],
        help='Different penalty values to test using cross validation.')

    # args for validating models:    
    validate_parent = argparse.ArgumentParser(add_help=False)
    validate_parent.add_argument("--validation_bam_file", help="hifi bam file", required=True)

    # args for validating models:    
    list_validate_parent = argparse.ArgumentParser(add_help=False)
    list_validate_parent.add_argument("--validation_list_file", help="file with list of true variants", required=True)


    # args for printing polished bam:    
    adjust_parent = argparse.ArgumentParser(add_help=False)
    adjust_parent.add_argument('--outbam', type=str,
        help="Bam file with adjusted base qualities.", required=True)
    adjust_parent.add_argument("--output_adjustments", type=argparse.FileType('w'))
    adjust_parent.add_argument("--cutoff", type=int,
        help="set adjusted basequalities lower than this cutoff to zero")

    # args for calling somatic variants:    
    call_parent = argparse.ArgumentParser(add_help=False)
    call_parent.add_argument('--outfile', type=argparse.FileType('w'), default=sys.stdout,
        help="output file")
    call_parent.add_argument('--method', type=str,
        choices=['LR', 'LR_with_MQ', 'poisson', 'BF'], default="LR",
        help="Method used to calculate variant quality scores")
    call_parent.add_argument('--cutoff', type=float, default=None, metavar='Q',
        help="Only print variants with quality above Q.")
    call_parent.add_argument('--prior_N', type=float, default=100,
        help="Weight (as sample size) of the kmer based prior on error rate.")
    call_parent.add_argument('--no_update',  action='store_true',
        help="Do not make bayesian update of error rate but use rate only estimated from kmers")
    call_parent.add_argument('--double_adjustment', type=float, default=0.5,
        help="How much lower should the error rate be if we see a matching overlap.")
    call_parent.add_argument('--min_BQ', type=int, default=1,
        help="Minimum base quality to considder")
    call_parent.add_argument('--min_MQ', type=int, default=50,
        help="Minimum base quality to considder")
    call_parent.add_argument('--filter_max_count', type=int, default=2,
        help='Maximum number of times an alternative read is allowed to be seen in filer_bam')
    call_parent.add_argument("--pop_vcf", type=str,
        help='Population vcf with AF field.')
    

    count_parser = subparsers.add_parser('count', 
        description='Count good and bad k-mers',
        help =  'Count good and bad k-mers',
        parents=[bam_parent, filter_parent, count_parent, read_filter_parent],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    train_parser = subparsers.add_parser('train', 
        description='First run "count" then train model to distinguish good and bad k-mers.', 
        help = 'First run "count" then train model to distinguish good and bad k-mers.',
        parents=[bam_parent, filter_parent, count_parent, train_parent, read_filter_parent])

    validate_parser = subparsers.add_parser('validate', 
        description = 'First run "count" and "train" then print validation data',
        help = 'First run "count" and "train" then print validation data',
        parents=[bam_parent, filter_parent, count_parent, train_parent, validate_parent, read_filter_parent])

    list_validate_parser = subparsers.add_parser('list_validate', 
        description = 'First run "count" and "train" then print validation data',
        help = 'First run "count" and "train" then print validation data',
        parents=[bam_parent, filter_parent, count_parent, train_parent, list_validate_parent, read_filter_parent])


    adjust_parser = subparsers.add_parser('adjust', 
        description = 'First run "count" and "train" then output bam with adjusted base qualities.', 
        help = 'First run "count" and "train" then output bam with adjusted base qualities.', 
        parents = [bam_parent, filter_parent, count_parent, train_parent, adjust_parent])

    call_parser = subparsers.add_parser('call', 
        description = 'First run "count" and "train" then call variants', 
        help = 'First run "count" and "train" then call variants', 
        parents = [bam_parent, count_parent, train_parent, call_parent, read_filter_parent])

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

    list_validate_only_parser = subparsers.add_parser('list_validate_only', 
        description = 'Print validation data.',
        help = 'Print validation data.', 
        parents = [bam_parent, filter_parent, list_validate_parent]) 
    list_validate_only_parser.add_argument("--input_file_kmerpapa", type=argparse.FileType('r'))

    adjust_only_parser = subparsers.add_parser('adjust_only', 
        description = 'Output bam with adjusted base qualities.', 
        help = 'Output bam with adjusted base qualities.',
        parents = [bam_parent, adjust_parent])
    adjust_only_parser.add_argument("--input_file_kmerpapa", type=argparse.FileType('r'))

    call_only_parser = subparsers.add_parser('call_only', 
        description = 'Call variants',
        help = 'Call variants',
        parents=[bam_parent, filter_parent, call_parent, read_filter_parent])
    call_only_parser.add_argument("--input_file_kmerpapa", type=argparse.FileType('r'))
    
    test_kmerpapa_parser = subparsers.add_parser('test_kmerpapa', 
        description = 'Apply a kmerpapa model to a set of kmer counts',
        help = 'Apply a kmerpapa model to a set of kmer counts')
    test_kmerpapa_parser.add_argument("--input_file_good", type=argparse.FileType('r'))
    test_kmerpapa_parser.add_argument("--input_file_bad", type=argparse.FileType('r'))
    test_kmerpapa_parser.add_argument("--input_file_kmerpapa", type=argparse.FileType('r'))
    test_kmerpapa_parser.add_argument('--correction_type', type=str, default = "bad_vs_no",
        help='should we compare bad variants to "good variants"(SNVs) or to "no variant" (homozygous ref sites)',
        choices=["bad_vs_good", "bad_vs_no"])
    test_kmerpapa_parser.add_argument('--same_good', action='store_true')
    return parser


def run_get_good_and_bad_w_filter(opts):
    if opts.verbosity > 0:
        eprint("Counting good and bad kmers")
    counter = \
        MutationCounterWFilter(
            opts.bam_file, 
            opts.twobit_file, 
            opts.filter_bam_file,
            radius = opts.radius,
            min_enddist = opts.min_enddist,
            max_mismatch = opts.max_mismatch)
    if opts.chrom is None:
        good_kmers, bad_kmers = \
            counter.count_mutations_all_chroms()
    else:
        good_kmers, bad_kmers = \
            counter.count_mutations(opts.chrom, opts.start, opts.end)
    print_good_and_bad(opts, good_kmers, bad_kmers)
    
    #change format of dicts to nested dicts
    BQs_good = set(x[0] for x in good_kmers)
    BQs_bad = set(x[0] for x in bad_kmers)
    good_kmers2 = {}
    bad_kmers2 = {}
    for BQ in BQs_good | BQs_bad:
        good_kmers2[BQ] = {}
        for mtype in ('A->C', 'A->G', 'A->T', 'C->A', 'C->G', 'C->T', 'C->C', 'A->A'):
            good_kmers2[BQ][mtype] = defaultdict(int)
        bad_kmers2[BQ] = {}
        for mtype in ('A->C', 'A->G', 'A->T', 'C->A', 'C->G', 'C->T'):
            bad_kmers2[BQ][mtype] = defaultdict(int)
    for tup, count in good_kmers.items():
        BQ, mtype, kmer = tup
        good_kmers2[BQ][mtype][kmer] = count
    for tup, count in bad_kmers.items():
        BQ, mtype, kmer = tup
        bad_kmers2[BQ][mtype][kmer] = count
    return good_kmers2, bad_kmers2


def run_get_kmerpapas(opts, good_kmers, bad_kmers):
    if opts.verbosity > 0:
        eprint("Training kmer pattern partitions")
    kmer_papas = {}

    for bqual in bad_kmers:
        radius = len(next(iter(good_kmers[bqual]["C->T"].keys())))//2
        kmer_papas[bqual] = {}
        eprint(f'Handling base_qual: {bqual}')
        for mtype in bad_kmers[bqual]:
            if mtype[0] == mtype[-1]:
                continue
            eprint(f'Handling mutation type: {mtype}')
            super_pattern = 'N'*radius + mtype[0] + 'N'*radius
            eprint(mtype, end=" ")
            if opts.correction_type == "bad_vs_good":
                if opts.same_good:
                    contextD = dict((x, (bad_kmers[bqual][mtype][x], good_kmers[37][mtype][x])) for x in matches(super_pattern))
                else:
                    contextD = dict((x, (bad_kmers[bqual][mtype][x], good_kmers[bqual][mtype][x])) for x in matches(super_pattern))
            elif opts.correction_type == "bad_vs_no":
                ref = mtype[0]
                alt = mtype[-1]
                notype = f'{ref}->{ref}'
                other_type1, other_type2 = [f'{ref}->{other}' for other in 'ACGT' if other != alt and other != ref]
                if opts.same_good:
                    contextD = dict((x, (bad_kmers[bqual][mtype][x], good_kmers[37][notype][x])) for x in matches(super_pattern))
                else:
                    contextD = dict((x, (bad_kmers[bqual][mtype][x], good_kmers[bqual][notype][x] + bad_kmers[bqual][other_type1][x] + bad_kmers[bqual][other_type2][x])) for x in matches(super_pattern))
            kpp = get_kmerpapa(super_pattern, contextD, opts)
            kmer_papas[bqual][mtype] = {}
            for pat in kpp:
                alpha, beta = kpp[pat]
                if not opts.output_file_kmerpapa is None:
                    print(bqual, mtype, pat, alpha, beta, file=opts.output_file_kmerpapa)
                if opts.verbosity > 1:
                    eprint(bqual, mtype, pat, alpha, beta, -10*log10(alpha/(alpha+beta)))
                for context in matches(pat):
                    kmer_papas[bqual][mtype][context] = (alpha, beta)
    
    if not opts.output_file_kmerpapa is None:
        opts.output_file_kmerpapa.close()
    
    return kmer_papas

def phred_scale_kmerpapas(kmer_papas):
    for bqual in kmer_papas:
        for mtype in kmer_papas[bqual]:
            for context in kmer_papas[bqual][mtype]:
                alpha, beta = kmer_papas[bqual][mtype][context]
                kmer_papas[bqual][mtype][context] = -10*log10(alpha/(alpha+beta))
    
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
    if opts.chrom is None:
        validator.call_all_chroms()
    else:
        validator.call_mutations(opts.chrom, opts.start, opts.end)

def run_list_validation(opts, kmer_papas):
    if opts.verbosity > 0:
        eprint("Printing validation data")
    validator = \
        ListMutationValidator(
            opts.bam_file, 
            opts.filter_bam_file, 
            opts.validation_list_file, 
            opts.twobit_file, 
            kmer_papas)
    if opts.chrom is None:
        validator.call_all_chroms()
    else:
        validator.call_mutations(opts.chrom, opts.start, opts.end)



def run_adjust(opts, kmer_papas):
    if opts.verbosity > 0:
        eprint("Adjusting base qualities")
    phred_scale_kmerpapas(kmer_papas)
    adjuster = \
        BaseAdjuster(
            opts.bam_file,
            opts.twobit_file, 
            kmer_papas,
            opts.outbam,
            opts.output_adjustments,
            opts.cutoff
        )
    
    n_corrections, n_corrected_reads, n_uncorrected_reads, n_filtered = \
        adjuster.call_mutations(opts.chrom, opts.start, opts.end)

    if opts.verbosity > 0:
        eprint(f'corrected {n_corrections} base qualities in {n_corrected_reads} reads')
        eprint(f'{n_uncorrected_reads} reads were written with no corrections')
        eprint(f'{n_filtered} reads were filtered')


def run_call(opts, kmer_papas):
    if opts.verbosity > 0:
        eprint("Calling somatic variants")
    caller = \
        SomaticMutationCaller(
            opts.bam_file,
            opts.filter_bam_file,
            opts.twobit_file, 
            kmer_papas,
            opts.outfile,
            opts.method,
            opts.cutoff,
            opts.prior_N,
            opts.no_update,
            opts.double_adjustment,
            opts.min_MQ,
            opts.min_BQ,
            opts.filter_max_count,
            opts.pop_vcf,
            opts.min_enddist,
            opts.max_mismatch
        )
    if opts.chrom is None:
        n_calls = caller.call_all_chroms()
    else:
        n_calls = caller.call_mutations(opts.chrom, opts.start, opts.end)

    if opts.verbosity > 0:
        eprint(f'Found {n_calls} variants with a p-value less than {opts.cutoff}.')
    

def get_phred(tup):
    alpha, beta = tup
    return -10*log10(alpha /(alpha+beta))
 
def run_test_kmerpapas(opts, kmer_papas, good_kmers, bad_kmers):
    for bqual in kmer_papas:
        for mtype in kmer_papas[bqual]:
            if mtype[0] == mtype[-1]:
                continue
            eprint(f'Handling base_qual: {bqual}, mutation type: {mtype}')
            for pat in kmer_papas[bqual][mtype]:
                if opts.correction_type == "bad_vs_good":
                    n_bad = sum(bad_kmers[bqual][mtype][x] for x in matches(pat))
                    if opts.same_good:
                        n_good = sum(good_kmers[37][mtype][x] for x in matches(pat))
                    else:
                        n_good = sum(good_kmers[bqual][mtype][x] for x in matches(pat))
                elif opts.correction_type == "bad_vs_no":
                    ref = mtype[0]
                    alt = mtype[-1]
                    notype = f'{ref}->{ref}'
                    other_type1, other_type2 = [f'{ref}->{other}' for other in 'ACGT' if other != alt and other != ref]
                    n_bad = sum(bad_kmers[bqual][mtype][x] for x in matches(pat))
                    if opts.same_good:
                        n_good = sum(good_kmers[37][notype][x] for x in matches(pat))
                    else:
                        n_good = sum(good_kmers[bqual][notype][x] + bad_kmers[bqual][other_type1][x] + bad_kmers[bqual][other_type2][x] for x in matches(pat))

                if n_bad > 0:
                    phred = -10*log10(n_bad / (n_bad + n_good))
                else:
                    phred = 0
                print(bqual, mtype, pat, get_phred(kmer_papas[bqual][mtype][pat]), n_bad, n_good, phred)


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

    if not opts.command in ['train_only', 'validate_only', 'list_validate_only', 'call_only', 'adjust_only', 'test_kmerpapa']:
        good_kmers, bad_kmers = run_get_good_and_bad_w_filter(opts)
    elif opts.command in ['train_only', 'test_kmerpapa']:
        good_kmers, bad_kmers = read_kmers(opts)

    if opts.command == 'count':
        return 0

    if not opts.command in ['validate_only', 'list_validate_only', 'call_only', 'adjust_only', 'test_kmerpapa']:
        kmer_papas = run_get_kmerpapas(opts, good_kmers, bad_kmers)
    elif opts.command in "test_kmerpapa":
        kmer_papas = read_kmer_papas_for_test(opts)
    else:
        eprint("Reading kmer pattern partitions")
        kmer_papas = read_kmer_papas(opts)

    if opts.command in ['train_only', 'train']:
        return 0

    if opts.command in ['validate', 'validate_only']:
        run_validation(opts, kmer_papas)
    elif opts.command in ['list_validate', 'list_validate_only']:
        run_list_validation(opts, kmer_papas)
    elif opts.command in ['call', 'call_only']:
        run_call(opts, kmer_papas)
    elif opts.command in ['adjust', 'adjust_only']:
        run_adjust(opts, kmer_papas)
    elif opts.command == 'test_kmerpapa':
        run_test_kmerpapas(opts, kmer_papas, good_kmers, bad_kmers)
    else:
        eprint("Unknown command: {opts.command}")
        return 1

    return 0


