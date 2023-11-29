import sys
from collections import Counter

ostrand = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def list_vals(x):
    return x[1:-1].split(',')

counter = Counter()

for line in sys.stdin:
    L = line.split('\t')
    
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
        newBQs = list_vals(D['newBQ'])
        assert(len(oldBQs)==len(newBQs))
        for i in range(len(oldBQs)):
            oldBQ = oldBQs[i]
            if ',' in oldBQ:
                max_oldBQ = max(int(x) for x in oldBQ[1:-1].split(','))
            counter[(muttype, oldBQs[i], max_oldBQ, newBQs[i], validation, D['n_mismatch'], D['N_A_37'])] += 1

columns = [
    'muttype',
    'oldBQ',
    'max_oldBQ',
    'newBQ',
    'validation',
    'n_mismatch',
    'NA37',
    'count',
]
print('\t'.join(columns))
for tup, count in counter.items():
    print('\t'.join(tup) + '\t' + str(count))
