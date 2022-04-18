import csv
import sys
from itertools import product, count
from operator import add
from functools import reduce, partial

mtypes = ('A->C', 'A->G', 'A->T', 'C->A', 'C->G', 'C->T')

complement = str.maketrans("ATCGN", "TAGCN")

def reverse_complement(s):
    return s.translate(complement)[::-1]

def regions(bed_file):
    with open(bed_file) as fp:
        reader = csv.reader(fp, delimiter="\t")
        for line in reader:
            yield line[0], int(line[1]), int(line[2])

def regions_if_sorted(bed_file):
    with open(bed_file) as fp:
        reader = csv.reader(fp, delimiter="\t")
        for line in reader:
            yield line[0], int(line[1]), int(line[2])

def generate_mutation_types(k):
    if k % 2 == 0:
        raise ValueError("k must be uneven")
    types = list()
    mut_pos = (k - 1) // 2
    for mer in map(
        partial(reduce, add),
        product(*(["ATGC"] * mut_pos + ["TC"] + ["ATGC"] * mut_pos)),
    ):
        alts = "ATG"
        if mer[mut_pos] == "T":
            alts = "AGC"
        for alt in alts:
            after = mer[:mut_pos] + alt + mer[mut_pos + 1 :]
            mut = ">".join([mer[mut_pos],alt])
            types.append("_".join([mut, mer]))
    return sorted(types)

mutation_types_3_mer = dict(zip(generate_mutation_types(3), count()))

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


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class Read:
    def __init__(self, pileup_read):
        pos = pileup_read.query_position
        # set attributes
        self.start = pileup_read.alignment.reference_start
        self.end = pileup_read.alignment.reference_end
        self.allel = pileup_read.alignment.query_sequence[pos]
        self.is_reverse = pileup_read.alignment.is_reverse
        self.base_qual = pileup_read.alignment.query_qualities[pos]
        self.query_name = pileup_read.alignment.query_name
        self.length = abs(pileup_read.alignment.template_length)
        #self.NH = pileup_read.alignment.get_tag("NH")

        # Process cigar stats
        cigar_stats = pileup_read.alignment.get_cigar_stats()[1]
        self.has_indel = sum(cigar_stats[1:4]) != 0
        self.has_clip = sum(cigar_stats[4:6]) != 0
        self.NM = cigar_stats[10]


class ReadPair:
    def __init__(self, read1, read2):
        self.length = read1.length
        self.allel = read1.allel
        self.start = read2.start if read1.is_reverse else read1.start
        self.min_base_qual = min(read1.base_qual, read2.base_qual)
        self.max_NM = max(read1.NM, read2.NM)
        self.has_indel = 1 if (read1.has_indel or read2.has_indel) else 0
        self.has_clip = 1 if (read1.has_clip or read2.has_clip) else 0


# def zip_pileups(p1, p2, p3):
#     while True:
#         try:
#             pileupcolumn1 = next(p1)
#             pileupcolumn2 = next(p2)
#             pileupcolumn3 = next(p3)
#             while(pileupcolumn1.reference_pos < max(pileupcolumn2.reference_pos, pileupcolumn3.reference_pos)):
#                 pileupcolumn1 = next(p1)
#             while(pileupcolumn2.reference_pos < max(pileupcolumn1.reference_pos, pileupcolumn3.reference_pos)):
#                 pileupcolumn2 = next(p2)
#             while(pileupcolumn3.reference_pos < max(pileupcolumn1.reference_pos, pileupcolumn2.reference_pos)):
#                 pileupcolumn3 = next(p3)
#             if pileupcolumn1.reference_pos == pileupcolumn2.reference_pos == pileupcolumn3.reference_pos:
#                 yield (pileupcolumn1, pileupcolumn2, pileupcolumn3)
#         except StopIteration:
#             break

#TODO: Check that this new general function match the commented one above
def zip_pileups(*pileups):
    while True:
        try:
            pileupcolumns = [None] * len(pileups)
            for i in range(len(pileups)):
                pileupcolumns[i]  = next(pileups[i])
            for i in range(len(pileups)):
                max_pos = max(x.reference_pos for x in pileupcolumns)
                while pileupcolumns[i].reference_pos < max_pos:
                    pileupcolumns[i] = next(pileups[i])
            if all(x.reference_pos == pileupcolumns[0].reference_pos for x in pileupcolumns):
                yield pileupcolumns
        except StopIteration:
            break
