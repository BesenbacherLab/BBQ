from betterbasequals.utils import eprint
from kmerpapa.algorithms import greedy_penalty_plus_pseudo
from kmerpapa.pattern_utils import get_M_U
from math import log10

def get_greedy_kmerpapa(super_pattern, contextD, opts):
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
    CV = greedy_penalty_plus_pseudo.BaysianOptimizationCV(super_pattern, contextD, opts.nfolds, opts.iterations, opts.seed)
    best_alpha, best_penalty, test_score = CV.get_best_a_c()

    if opts.verbosity > 0:
        eprint(f'best_penalty = {best_penalty}, best_alpha={best_alpha}')
    best_beta = (best_alpha*(1.0-my))/my
    best_score, M, U, names = greedy_penalty_plus_pseudo.greedy_partition(super_pattern, contextD, best_alpha, best_beta, best_penalty, {})
    counts = []
    for pat in names:
        counts.append(get_M_U(pat, contextD))
    D = {}
    for i in range(len(names)):
        pat = names[i]
        M, U = counts[i]
        p = (M + best_alpha)/(M + U + best_alpha + best_beta)
        D[pat] = -10*log10(p)
    return D

def get_optimal_kmerpapa(super_pattern, contextD, opts):
    pseudo_counts = [1,2,5,10,20,30,50,100,500,1000]
    penalty_values = range(1,15)
    
