import sys

ostrand = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

for line in sys.stdin:
    L = line.split('\t')
    
    fromA = L[3]
    toA = L[4]

    if fromA in ['T','G']:
        muttype=f'{ostrand[fromA]}->{ostrand[toA]}'
    else:
        muttype=f'{fromA}->{toA}'

    if L[6] == 'PASS' or L[6] == 'lowBQ':
        D = dict(x.split('=') for x in  L[7].split(';'))
        print(line.strip())
        print(muttype)
        print(D)

