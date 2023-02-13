from betterbasequals.utils import eprint
from kmerpapa.algorithms import greedy_penalty_plus_pseudo
from kmerpapa.algorithms import bottum_up_array_penalty_plus_pseudo_CV
from kmerpapa.algorithms import bottum_up_array_w_numba
from kmerpapa.pattern_utils import get_M_U, generality
from math import log10
import argparse

def get_kmerpapa(super_pattern, contextD, opts):
    n_bad = 0
    n_good = 0
    for x,y in contextD.values():
        n_bad += x
        n_good += y
    my=n_bad/(n_good+n_bad)
    if opts.verbosity > 0:
        eprint(f'n_bad = {n_bad}, n_good = {n_good}, mean_error_rate = {my}')
    if n_bad<10:
        assert n_good > 50, f'too few counts'
        if opts.verbosity > 0:
            eprint('Skipping kmerpapa because of low number of bad kmers')
        return {super_pattern: -10*log10((0.5+n_bad)/(0.5+n_good +n_bad))}
    if n_good<10:
        assert n_bad > 50, f'too few counts'
        if opts.verbosity > 0:
            eprint('skipping kmerpapa because of low number of good kmers')
        return {super_pattern: -10*log10((n_bad)/(0.5+n_good +n_bad))}
    if opts.kmerpapa_method == 'greedy':
        CV = greedy_penalty_plus_pseudo.BaysianOptimizationCV(super_pattern, contextD, opts.nfolds, opts.iterations, opts.seed)
        best_alpha, best_penalty, test_score = CV.get_best_a_c()
        eprint(test_score)
        best_beta = (best_alpha*(1.0-my))/my
        best_score, M, U, names = greedy_penalty_plus_pseudo.greedy_partition(super_pattern, contextD, best_alpha, best_beta, best_penalty, opts)
    elif opts.kmerpapa_method == 'optimal':
        args = argparse.Namespace()
        args.nfolds = opts.nfolds
        args.iterations = opts.iterations
        args.verbosity = opts.verbosity
        args.CVfile = None
        args.seed = opts.seed
        pseudo_counts = opts.pseudo_counts
        penalty_values = opts.penalty_values
        best_alpha, best_penalty, test_score = bottum_up_array_penalty_plus_pseudo_CV.pattern_partition_bottom_up(super_pattern, contextD, pseudo_counts, args, n_bad, n_good, penalty_values)
        best_beta = (best_alpha*(1.0-my))/my
        best_score, M, U, names = bottum_up_array_w_numba.pattern_partition_bottom_up(super_pattern, contextD, best_alpha, best_beta, best_penalty, opts, n_bad, n_good, index_mut=0)
    if opts.verbosity > 0:
        eprint(f'best_penalty = {best_penalty}, best_alpha={best_alpha}, n_patterns={len(names)}')
    counts = []
    for pat in names:
        counts.append(get_M_U(pat, contextD))
    D = {}
    x = 0
    for i in range(len(names)):
        pat = names[i]
        M, U = counts[i]
        p = (M + best_alpha)/(M + U + best_alpha + best_beta)
        eprint(pat, p)
        D[pat] = -10*log10(p)
        x += generality(pat)
    assert x == generality(super_pattern)
    return D
