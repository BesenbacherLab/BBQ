"""Module that contains the command line application."""

import argparse
from betterbasequals.get_good_bad_kmers import MutationCounterWFilter
from betterbasequals.model_validators import MutationValidator, ListMutationValidator
from betterbasequals.bam_adjusters import BaseAdjuster
from betterbasequals.somatic_callers import SomaticMutationCaller
from betterbasequals.utils import *
from betterbasequals import __version__
from betterbasequals.kmerpapa_utils import get_kmerpapa
from betterbasequals.filter_calls import BBQFilter
from math import sqrt

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
    read_filter_parent.add_argument('--min_enddist', type=int, default=5, metavar="M",
        help="Ignore bases in the first M or last M positions in the read")
    read_filter_parent.add_argument('--max_mismatch', type=int, default=1, metavar="M",
        help="Ignore alt reads if the read has more than M mismatches to the reference")
    read_filter_parent.add_argument('--min_MQ', type=int, default=40,
        help="Minimum base quality to considder")
    read_filter_parent.add_argument('--exclude_bed', type=str,
        help='Bed file with regions to exclude. Should be sorted by chrom (alphabetically) then pos (numrically).')

    # args for counting kmers:
    count_parent = argparse.ArgumentParser(add_help=False)
    count_parent.add_argument("--output_file_kmers", type=argparse.FileType('w'))
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
    train_parent.add_argument("--output_file_EQ", type=argparse.FileType('w'))
    train_parent.add_argument("--CVfile_prefix", type=str)
    train_parent.add_argument('-N', '--nfolds', type=int, metavar='N', default=2,
        help='Number of folds to use when fitting hyperparameters in kmerpapa')
    train_parent.add_argument('-i', '--iterations', type=int, default=1, metavar='i',
        help='Repeat cross validation i times when fitting hyperparameters in kmerpapa')
    train_parent.add_argument('--seed', type=int,
        help='seed for numpy.random')
    train_parent.add_argument('--same_good', action='store_true')
    train_parent.add_argument('--subtract', action='store_true')
    train_parent.add_argument('--EQ_pat', action='store_true')
    train_parent.add_argument('--EQ_pat_alpha', type=float, default=2)
    train_parent.add_argument('--min1EQ', action='store_true')
    train_parent.add_argument('--estimated', type=str, default = 'single',
                              choices = ['single', 'double'])
    train_parent.add_argument('--mean_type', type=str, default = 'geometric',
                              choices = ['geometric', 'harmonic'])
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
        choices=['LR', 'LR_with_MQ', 'maxLR_with_MQ', 'BF', 'BF_with_MQ', 'BF_with_MQ_and_Prior'], default="LR",
        help="Method used to calculate variant quality scores")
    call_parent.add_argument('--cutoff', type=float, default=None, metavar='Q',
        help="Only print variants with quality above Q.")
    call_parent.add_argument('--prior_N', type=float, default=500,
        help="Weight (as sample size) of the kmer based prior on error rate.")
    call_parent.add_argument('--no_update',  action='store_true',
        help="Do not make bayesian update of error rate but use rate only estimated from kmers")
    call_parent.add_argument('--min_BQ', type=int, default=1,
        help="Minimum base quality to considder")
    call_parent.add_argument('--filter_max_count', type=int, default=2,
        help='Maximum number of times an alternative read is allowed to be seen in filer_bam')
    call_parent.add_argument("--pop_vcf", type=str,
        help='Population vcf with AF field.')
    #call_parent.add_argument('--mean_type', type=str,
    #    choices=['arithmetric', 'geometric'], default="geometric",
    #    help="How to calculate the mean error probabilty of two matching overlapping alleles")
    call_parent.add_argument('--max_NM_diff', type=int, default=1,
        help="Maximum allowed difference in median number of mismatches for reads with alt alleles compared to reads with ref alleles")

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
    train_only_parser.add_argument("--input_file_kmers", type=argparse.FileType('r'))
    #train_only_parser.add_argument("--input_file_bad", type=argparse.FileType('r'))

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
    test_kmerpapa_parser.add_argument("--input_file_kmers", type=argparse.FileType('r'))
    #test_kmerpapa_parser.add_argument("--input_file_bad", type=argparse.FileType('r'))
    test_kmerpapa_parser.add_argument("--input_file_kmerpapa", type=argparse.FileType('r'))
    test_kmerpapa_parser.add_argument('--correction_type', type=str, default = "bad_vs_no",
        help='should we compare bad variants to "good variants"(SNVs) or to "no variant" (homozygous ref sites)',
        choices=["bad_vs_good", "bad_vs_no"])
    test_kmerpapa_parser.add_argument('--same_good', action='store_true')

    filter_calls_parser = subparsers.add_parser('filter_calls', 
        description = 'Filter PASS calls based on coverage quantiles',
        help = 'Filter PASS calls based on coverage quantiles')
    filter_calls_parser.add_argument('--vcf_file', required=True,
        help='input VCF file')
    filter_calls_parser.add_argument('--outfile',  type=argparse.FileType('w'), default=sys.stdout,
        help='output VCF file') 
    filter_calls_parser.add_argument('--lower_q', type=float, default=0.05,
        help='Lower quantile, PASS calls with coverage below this value will be filtered out') 
    filter_calls_parser.add_argument('--upper_q', type=float, default=0.95,
        help='Upper quantile, PASS calls with coverage above this value will be filtered out') 
    return parser


def run_get_good_and_bad_w_filter(opts):
    if opts.verbosity > 0:
        eprint("Counting good and bad kmers")
    counter = \
        MutationCounterWFilter(
            opts.bam_file, 
            opts.twobit_file, 
            opts.filter_bam_file,
            opts.min_MQ,
            radius = opts.radius,
            min_enddist = opts.min_enddist,
            max_mismatch = opts.max_mismatch,
            bed_file = opts.exclude_bed)
    
    if opts.chrom is None:
        event_kmers = counter.count_mutations_all_chroms()
    else:
        event_kmers = counter.count_mutations(opts.chrom, opts.start, opts.end)
    print_kmer_counts(opts, event_kmers)
    
    return event_kmers
    
    #change format of dicts to nested dicts
    # BQs = set(x[1] for x in event_kmers)
    # good_kmers = {}
    # bad_kmers = {}
    # for BQ in BQs:
    #     good_kmers[BQ] = {}
    #     for mtype in ('A->C', 'A->G', 'A->T', 'C->A', 'C->G', 'C->T', 'C->C', 'A->A'):
    #         good_kmers[BQ][mtype] = defaultdict(int)
    #     bad_kmers[BQ] = {}
    #     for mtype in ('A->C', 'A->G', 'A->T', 'C->A', 'C->G', 'C->T'):
    #         bad_kmers[BQ][mtype] = defaultdict(int)

    # for tup, count in event_kmers.items():
    #     event_type, BQ, mtype, kmer = tup
    #     if event_type =='good':
    #         good_kmers[BQ][mtype][kmer] = count
    #     elif event_type == 'bad':
    #         bad_kmers[BQ][mtype][kmer] = count
    # return good_kmers, bad_kmers


def run_get_kmerpapas(opts, event_kmers):
    if opts.verbosity > 0:
        eprint("Training kmer pattern partitions")
    kmer_papas = {}
    kmer_patterns = {}
    BQs = list(set(x[1] for x in event_kmers if x[0]in ['good', 'bad']))
    BQs.sort()
    BQs = BQs[::-1]
    radius = len(next(iter(event_kmers.keys()))[3])//2
    for bqual in BQs:
        BQ_pair = f'({bqual},{bqual})'
        kmer_papas[BQ_pair] = {}
        kmer_papas[bqual] = {}
        kmer_patterns[bqual] = {}
        for mtype in ('A->C', 'A->G', 'A->T', 'C->A', 'C->G', 'C->T'):
            contextD = {}
            ref = mtype[0]
            alt = mtype[-1]
            super_pattern = 'N'*radius + ref + 'N'*radius
            notype = f'{ref}->{ref}'
            other_type1, other_type2 = [f'{ref}->{other}' for other in 'ACGT' if other != alt and other != ref]
            eprint(f'Handling base_qual: {BQ_pair}, muttype: {mtype}')
            for kmer in matches(super_pattern):
                n_errors = event_kmers[('bad', bqual, mtype, kmer)]
                n_not_error = event_kmers[('good', bqual, notype, kmer)]
                if not opts.same_good:
                    n_not_error += event_kmers[('bad', bqual, other_type1, kmer)]
                    n_not_error += event_kmers[('bad', bqual, other_type2, kmer)]
                contextD[kmer] = (n_errors, n_not_error)

            kpp = get_kmerpapa(super_pattern, contextD, f'{mtype}_{bqual}', opts)
            if opts.estimated == "single":
                estBQ = bqual
            elif opts.estimated == "double":
                estBQ = BQ_pair

            kmer_papas[bqual][mtype] = {}    
            kmer_papas[BQ_pair][mtype] = {}
            
            kmer_patterns[bqual][mtype] = kpp#list(kpp.keys())
            for pat in kpp:
                alpha, beta = kpp[pat]
                p_error = min(0.25, alpha / (alpha + beta))
                if opts.verbosity > 1:
                    eprint(estBQ, mtype, pat, p_error, -10*log10(p_error))
                for kmer in matches(pat):
                    kmer_papas[estBQ][mtype][kmer] = p_error

    # Estimate EQ correction for bases in overlaps
    single_EQ = {}
    double_EQ = {}
    avg_single_EQ = {}
    avg_double_EQ = {}
    for BQ in BQs:
        single_EQ[BQ] = {}
        double_EQ[BQ] = {}
        avg_single_EQ[BQ] = {}
        avg_double_EQ[BQ] = {}
        for mtype in ('A->C', 'A->G', 'A->T', 'C->A', 'C->G', 'C->T'):
            ref = mtype[0]
            alt = mtype[-1]
            super_pattern = 'N'*radius + ref + 'N'*radius
            notype = f'{ref}->{ref}'
            other_type1, other_type2 = [f'{ref}->{other}' for other in 'ACGT' if other != alt and other != ref]
            BQ_pair = f'({BQ},{BQ})'
            single_mut = 0
            single_nomut = 0
            double_mut = 0
            double_nomut = 0
            if opts.estimated == 'single':
                BQest = BQ
            elif opts.estimated == 'double':
                BQest = BQ_pair
            
            if opts.EQ_pat:
                total_single_mut = 0
                total_single_nomut = 0
                total_double_mut = 0
                total_double_nomut = 0
                single_EQ[BQ][mtype] = {}
                double_EQ[BQ][mtype] = {}
                avg_single_EQ[BQ][mtype] = {}
                avg_double_EQ[BQ][mtype] = {}
                for pat in kmer_patterns[BQ][mtype]:
                    single_mut = 0
                    single_nomut = 0
                    double_mut = 0
                    double_nomut = 0
                    for kmer in matches(pat):
                        single_mut += event_kmers[('singleton', BQ, mtype, kmer)]
                        single_nomut += event_kmers[('singleton', BQ, notype, kmer)]
                        double_mut += event_kmers[('good_tuple', BQ_pair, mtype, kmer)]
                        double_nomut += event_kmers[('good_tuple', BQ_pair, notype, kmer)]
                        #TODO: add other muttypes from good_tuple and single to nomut?
                        double_nomut += event_kmers[('bad_tuple', BQ_pair, mtype, kmer)]
                        double_nomut += event_kmers[('bad_tuple', BQ_pair, other_type1, kmer)]
                        double_nomut += event_kmers[('bad_tuple', BQ_pair, other_type2, kmer)]
                    total_single_mut += single_mut
                    total_single_nomut += single_nomut
                    total_double_mut += double_mut
                    total_double_nomut += double_nomut
                    single_EQ[BQ][mtype][pat] = (single_mut, single_nomut)#(single_mut+0.1) / (single_mut + single_nomut + 0.2)
                    double_EQ[BQ][mtype][pat] = (double_mut, double_nomut)#(double_mut+0.1) / (double_mut + double_nomut + 0.2)
                avg_single_EQ[BQ][mtype] = (total_single_mut + 0.1) / (total_single_mut + total_single_nomut + 0.2)
                avg_double_EQ[BQ][mtype] = (total_double_mut + 0.1) / (total_double_mut + total_double_nomut + 0.2)

            else:
                for kmer in kmer_papas[BQest][mtype]:
                    single_mut += event_kmers[('singleton', BQ, mtype, kmer)]
                    single_nomut += event_kmers[('singleton', BQ, notype, kmer)]
                    double_mut += event_kmers[('good_tuple', BQ_pair, mtype, kmer)]
                    double_nomut += event_kmers[('good_tuple', BQ_pair, notype, kmer)]
                    #TODO: add other muttypes from good_tuple and single to nomut?
                    double_nomut += event_kmers[('bad_tuple', BQ_pair, mtype, kmer)]
                    double_nomut += event_kmers[('bad_tuple', BQ_pair, other_type1, kmer)]
                    double_nomut += event_kmers[('bad_tuple', BQ_pair, other_type2, kmer)]
                single_EQ[BQ][mtype] = (single_mut+0.1) / (single_mut + single_nomut + 0.2)
                double_EQ[BQ][mtype] = (double_mut+0.1) / (double_mut + double_nomut + 0.2)
           
           # rel_EQ = single_EQ / double_EQ

            #if not opts.output_file_EQ is None:
            #    print(BQ, mtype, single_mut, single_nomut, single_EQ, double_mut, double_nomut, double_EQ, rel_EQ, file = opts.output_file_EQ)

    #if not opts.output_file_EQ is None:
    #    opts.output_file_EQ.close()        

    if not opts.output_file_EQ is None:
        if opts.EQ_pat:
            print('BQ mutationtype pattern kmerpapa_alpha kmepapa_beta kmerpapa_rate single_mut single_nomut double_mut double_nomut single_EQ double_EQ subtract relative_EQ', file = opts.output_file_EQ)



    max_BQ = max(BQs)
    alpha = opts.EQ_pat_alpha #5 # pseudo count
    # apply EQ correction
    for BQ in BQs:
        for mtype in ('A->C', 'A->G', 'A->T', 'C->A', 'C->G', 'C->T'):
            BQ_pair = f'({BQ},{BQ})'
            if (opts.subtract and
                (BQ < max_BQ) and 
                (single_EQ[BQ][mtype] > double_EQ[max_BQ][mtype] + 1e-8) and 
                (double_EQ[BQ][mtype] > double_EQ[max_BQ][mtype] + 1e-8)):
                subtract = double_EQ[max_BQ][mtype]
            else:
                subtract = 0.0

            if not opts.EQ_pat:
                rel_EQ = (single_EQ[BQ][mtype] - subtract) / (double_EQ[BQ][mtype] - subtract)
                if not opts.output_file_EQ is None:
                    print(BQ, mtype, single_EQ[BQ][mtype], double_EQ[BQ][mtype], subtract, rel_EQ, file = opts.output_file_EQ)

            if opts.estimated == "single":
                kmer_papas[BQ_pair][mtype] = {}
                if opts.EQ_pat:
                    for pat in single_EQ[BQ][mtype]:
                        single_mut, single_nomut = single_EQ[BQ][mtype][pat]
                        double_mut, double_nomut =  double_EQ[BQ][mtype][pat]    
                        beta_single = (alpha*(1.0-avg_single_EQ[BQ][mtype]))/avg_single_EQ[BQ][mtype]
                        beta_double = (alpha*(1.0-avg_double_EQ[BQ][mtype]))/avg_double_EQ[BQ][mtype]
                        s_EQ = (single_mut + alpha) / (single_mut + single_nomut + alpha +beta_single)
                        d_EQ = (double_mut + alpha) / (double_mut + double_nomut + alpha +beta_double)
                        rel_EQ = s_EQ / d_EQ
                        if opts.min1EQ:
                            rel_EQ = max(1.0, rel_EQ)
                        if not opts.output_file_EQ is None:
                            kmp_alpha, kmp_beta = kmer_patterns[BQ][mtype][pat]
                            print(BQ, mtype, pat, kmp_alpha, kmp_beta, kmp_alpha/(kmp_alpha+kmp_beta), single_mut, single_nomut, double_mut, double_nomut, s_EQ, d_EQ, subtract, rel_EQ, file = opts.output_file_EQ)
                        for kmer in matches(pat):
                            kmer_papas[BQ_pair][mtype][kmer] = min(kmer_papas[BQ][mtype][kmer] / rel_EQ, 0.25)
                
                else:
                    for kmer in kmer_papas[BQ][mtype]:
                        kmer_papas[BQ_pair][mtype][kmer] = min(kmer_papas[BQ][mtype][kmer] / rel_EQ, 0.25)
                

            elif opts.estimated == "double":
                kmer_papas[BQ][mtype] = {}
                for kmer in kmer_papas[BQ_pair][mtype]:
                    kmer_papas[BQ][mtype][kmer] = min(kmer_papas[BQ_pair][mtype][kmer] * rel_EQ, 0.25)
    
    if not opts.output_file_EQ is None:
        opts.output_file_EQ.close() 

    # Estimate error rates for pairs of different BQs:
    for i in range(len(BQs)):
        for j in range(i+1, len(BQs)):
            BQ1 = BQs[i]
            BQ2 = BQs[j]
            BQ_pair = f'({BQ1},{BQ2})'
            BQ1_pair = f'({BQ1},{BQ1})'
            BQ2_pair = f'({BQ2},{BQ2})'
            kmer_papas[BQ_pair] = {}
            for mtype in ('A->C', 'A->G', 'A->T', 'C->A', 'C->G', 'C->T'):
                kmer_papas[BQ_pair][mtype] = {}
                for kmer in kmer_papas[BQ1_pair][mtype]:
                    if opts.mean_type == 'geometric':
                        kmer_papas[BQ_pair][mtype][kmer] = \
                            sqrt(kmer_papas[BQ1_pair][mtype][kmer]*kmer_papas[BQ2_pair][mtype][kmer])
                    elif opts.mean_type == 'harmonic':
                        kmer_papas[BQ_pair][mtype][kmer] = \
                            1.0/(((1.0/kmer_papas[BQ1_pair][mtype][kmer]) + (1.0/kmer_papas[BQ2_pair][mtype][kmer]))/2)

                    #TODO: Should we try harmonic mean?

    if not opts.output_file_kmerpapa is None:
        for BQ in kmer_papas:
            for mtype in kmer_papas[BQ]:
                for kmer in kmer_papas[BQ][mtype]:
                    print(BQ, mtype, kmer, kmer_papas[BQ][mtype][kmer], file=opts.output_file_kmerpapa)                   
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
    if opts.prior_N < 0:
        opts.no_update = True
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
            opts.min_MQ,
            opts.min_BQ,
            opts.filter_max_count,
            opts.pop_vcf,
            opts.min_enddist,
            opts.max_mismatch,
            opts.max_NM_diff,
            bed_file = opts.exclude_bed,
        )
    if opts.chrom is None:
        n_calls = caller.call_all_chroms()
    else:
        n_calls = caller.call_mutations(opts.chrom, opts.start, opts.end)

    if opts.verbosity > 0:
        if opts.cutoff is None:
            eprint(f'Found {n_calls} possible variants.')
        else:
            eprint(f'Found {n_calls} variants with a quality above {opts.cutoff}.')
    

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


def run_calls_filter(opts):
    if opts.verbosity > 0:
        eprint(f'Filtering PASS calls based on {opts.lower_q*100}% and {opts.upper_q*100}% coverage quantiles.')

    filter = \
        BBQFilter(
            opts.vcf_file, 
            opts.outfile, 
            opts.lower_q,
            opts.upper_q
            )
    n_filtered = filter.filter_BBQ()
    
    if opts.verbosity > 0:
        eprint(f"Removed {n_filtered} PASS calls.")


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

    if not opts.command in ['train_only', 'validate_only', 'list_validate_only', 'call_only', 'adjust_only', 'test_kmerpapa', 'filter_calls']:
        event_kmers = run_get_good_and_bad_w_filter(opts)
    elif opts.command in ['train_only', 'test_kmerpapa']:
        event_kmers = read_kmers(opts)

    if opts.command == 'count':
        return 0

    if not opts.command in ['validate_only', 'list_validate_only', 'call_only', 'adjust_only', 'test_kmerpapa', 'filter_calls']:
        kmer_papas = run_get_kmerpapas(opts, event_kmers)
    elif opts.command in "test_kmerpapa":
        kmer_papas = read_kmer_papas_for_test(opts)
    elif opts.command in "filter_calls":
        pass
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
    #elif opts.command == 'test_kmerpapa':
    #    run_test_kmerpapas(opts, kmer_papas, event_kmers)
    elif opts.command == 'filter_calls':
        run_calls_filter(opts)
    else:
        eprint("Unknown command: {opts.command}")
        return 1

    return 0


