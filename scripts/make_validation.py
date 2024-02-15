import sys
from collections import Counter
import argparse

complement = str.maketrans("ATCGN", "TAGCN")

def reverse_complement(s):
    return s.translate(complement)[::-1]

ostrand = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

code = {'A':['A'],
        'C':['C'],
        'G':['G'],
        'T':['T'],
        'R':['A', 'G'],
        'Y':['C', 'T'],
        'S':['G', 'C'],
        'W':['A', 'T'],
        'K':['G', 'T'],
        'M':['A', 'C'],
        'B':['C', 'G', 'T'],
        'D':['A', 'G', 'T'],
        'H':['A', 'C', 'T'],
        'V':['A', 'C', 'G'],
        'N':['A', 'C', 'G', 'T']}

def matches(pattern):
    if len(pattern) == 0:
        yield ''
    else:
        for y in matches(pattern[1:]):
            for x in code[pattern[0]]:
                yield x+y


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
        L = line.split()
        BQ = L[0]
        muttype = L[1]
        pattern = L[2]
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
    

    pos_type = f'{fromA}->{toA}'
    neg_type = f'{ostrand[fromA]}->{ostrand[toA]}'
    if fromA in ['T','G']:
        sym_muttype=f'{ostrand[fromA]}->{ostrand[toA]}'
    else:
        sym_muttype=f'{fromA}->{toA}'

    if L[6] == 'PASS' or L[6] == 'lowBQ':
        D = dict(x.split('=') for x in  L[7].strip().split(';'))
        oldBQs = list_vals(D['oldBQ'])
        newBQs = [str(int(float(x))) for x in list_vals(D['newBQ'])]
        strands = list_vals(D['strand'])
        kmer = D['kmer']
        
        assert(len(oldBQs)==len(newBQs))
        for i in range(len(oldBQs)):
            oldBQ = oldBQs[i]
            strand = strands[i]
            if strand == '0':
                nonsym_muttype = pos_type
            elif strand == '1':
                nonsym_muttype = neg_type
            else:
                assert(False)
            
            if '/' in oldBQ:
                max_oldBQ = str(max(int(x) for x in oldBQ.split('/')))
            else:
                max_oldBQ = oldBQ
            if kmer2pattern is None:
                counter[(muttype, nonsym_muttype, oldBQs[i], max_oldBQ, newBQs[i], validation, D['n_mismatch'], D['N_A_37'])] += 1
            else:
                if kmer not in kmer2pattern[max_oldBQ][muttype]:
                    kmer = reverse_complement(kmer)
                pattern = kmer2pattern[max_oldBQ][muttype][kmer]
                counter[(muttype, nonsym_muttype, oldBQs[i], max_oldBQ, newBQs[i], validation, D['n_mismatch'], D['N_A_37'], pattern)] += 1

columns = [
    'sym_muttype',
    'nonsym_muttype',
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
