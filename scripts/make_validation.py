import sys
from collections import Counter
from betterbasequals.utils import matches, reverse_complement
import argparse

ostrand = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def list_vals(x):
    return x[1:-1].split(',')

counter = Counter()

parser = argparse.ArgumentParser(description='''

''')
parser.add_argument('--EQ_file', type=argparse.FileType('r'))
args = parser.parse_args()

kmer2pattern = None
if not args.EQ_file is None:
    kmer2pattern = {}
    pattern2correction = {}
    for line in args.EQ_file:
        BQ, muttype, pattern, single_rate, double_rat, subtract, correction = line.split()
        if BQ not in kmer2pattern:
            kmer2pattern[BQ] = {}
        if muttype not in kmer2pattern[BQ]:
            kmer2pattern[BQ][muttype] = {}
        for kmer in matches(pattern):
            kmer2pattern[BQ][muttype][kmer] = pattern


for line in sys.stdin:
    L = line.strip().split('\t')
    
    fromA = L[3]
    toA = L[4]
    validation = L[-1]

    if fromA in ['T','G']:
        muttype=f'{ostrand[fromA]}->{ostrand[toA]}'
    else:
        muttype=f'{fromA}->{toA}'

    if L[6] == 'PASS' or L[6] == 'lowBQ':
        D = dict(x.split('=') for x in  L[7].strip().split(';'))
        oldBQs = list_vals(D['oldBQ'])
        newBQs = [str(int(float(x))) for x in list_vals(D['newBQ'])]
        kmer = D['kmer']
        
        assert(len(oldBQs)==len(newBQs))
        for i in range(len(oldBQs)):
            oldBQ = oldBQs[i]
            if '/' in oldBQ:
                max_oldBQ = str(max(int(x) for x in oldBQ.split('/')))
            else:
                max_oldBQ = oldBQ
            if kmer2pattern is None:
                counter[(muttype, oldBQs[i], max_oldBQ, newBQs[i], validation, D['n_mismatch'], D['N_A_37'])] += 1
            else:
                if kmer not in kmer2pattern[max_oldBQ][muttype]:
                    kmer = reverse_complement(kmer)
                pattern = kmer2pattern[max_oldBQ][muttype][kmer]
                counter[(muttype, oldBQs[i], max_oldBQ, newBQs[i], validation, D['n_mismatch'], D['N_A_37'], pattern)] += 1

columns = [
    'muttype',
    'oldBQ',
    'max_oldBQ',
    'newBQ',
    'validation',
    'n_mismatch',
    'NA37',
]

if not kmer2pattern is None:
    columns.append('pattern')
columns.append('count')

print('\t'.join(columns))
for tup, count in counter.items():
    print('\t'.join(tup) + '\t' + str(count))
